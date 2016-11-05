#include <cassert>
#include <mutex>
#include <limits>
#include <vector>
#include <string>
#include <thread>
#include <fstream>
#include <iostream>
#include <exception>
#include <algorithm>
#include <unordered_map>
#include <condition_variable>

using namespace std;

#include <MpiSupport.h>
#include <Vector2d.h>
#include <DissimilarityMatrix.h>
#include <PartitioningAroundMedoids.h>

typedef float DistanceType;
typedef CVector2d<DistanceType> CVector;
typedef CDissimilarityMatrix<DistanceType> DissimilarityMatrixType;
typedef CPartitioningAroundMedois<DissimilarityMatrixType> PamType;

///////////////////////////////////////////////////////////////////////////////

class CBarrier {
	CBarrier( const CBarrier& ) = delete;
	CBarrier& operator=( const CBarrier& ) = delete;

public:
	explicit CBarrier( size_t _numberOfThreads ) :
		numberOfThreads( _numberOfThreads ),
		counter( 0 ),
		isCountingDown( false )
	{
		if( numberOfThreads == 0 ) {
			throw invalid_argument( "CBarrier: number of threads must be positive" );
		}
	}

	void Sync()
	{
		unique_lock<mutex> lock{ m };

		if( isCountingDown ) {
			counter--;
			if( counter == 0 ) {
				isCountingDown = false;
				cv.notify_all();
			} else {
				cv.wait( lock, [this]{ return !isCountingDown; } );
			}
		} else {
			counter++;
			if( counter == numberOfThreads ) {
				isCountingDown = true;
				cv.notify_all();
			} else {
				cv.wait( lock, [this]{ return isCountingDown; } );
			}
		}
	}

private:
	mutex m;
	condition_variable cv;

	size_t counter;
	const size_t numberOfThreads;
	bool isCountingDown;
};

///////////////////////////////////////////////////////////////////////////////

static_assert( sizeof( size_t ) <= sizeof( uint32_t ),
	"invalid: sizeof( size_t ) <= sizeof( uint32_t )" );

struct CObjectMedoidDistance {
	uint32_t Object;
	uint32_t Medoid;
	DistanceType Distance;

	CObjectMedoidDistance() :
		Object( 0 ),
		Medoid( 0 ),
		Distance( 0 )
	{
	}

	void Min( const CObjectMedoidDistance& another );
	void AllReduce();

private:
	static MPI_Datatype datatype();
	static MPI_Op op();
	static void MPIAPI objectMedoidDistanceMin(
		CObjectMedoidDistance* in, CObjectMedoidDistance* inout,
		int* length, MPI_Datatype* /*type*/ );
};

void CObjectMedoidDistance::Min( const CObjectMedoidDistance& another )
{
	if( another.Distance < Distance ) {
		*this = another;
	}
}

void CObjectMedoidDistance::AllReduce()
{
	MpiCheck( MPI_Allreduce( MPI_IN_PLACE, this, 1 /* count */, datatype(), op(), MPI_COMM_WORLD ),
		"MPI_Allreduce for CObjectMedoidDistance" );
}

MPI_Datatype CObjectMedoidDistance::datatype()
{
	static MPI_Datatype type = MPI_DATATYPE_NULL;
	if( type == MPI_DATATYPE_NULL ) {
		const int count = 3;
		int blocklengths[count] = { 1, 1, 1 };
		MPI_Datatype types[count] = { MPI_UINT32_T, MPI_UINT32_T, MPI_FLOAT };
		MPI_Aint offsets[count] = {
			offsetof( CObjectMedoidDistance, Object ),
			offsetof( CObjectMedoidDistance, Medoid ),
			offsetof( CObjectMedoidDistance, Distance ) };
		MpiCheck( MPI_Type_create_struct( count, blocklengths, offsets, types, &type ),
			"MPI_Type_create_struct for CObjectMedoidDistance" );
		MpiCheck( MPI_Type_commit( &type ), "MPI_Type_commit for CObjectMedoidDistance" );
	}
	return type;
}

MPI_Op CObjectMedoidDistance::op()
{
	static MPI_Op op = MPI_OP_NULL;
	if( op == MPI_OP_NULL ) {
		MpiCheck( MPI_Op_create( (MPI_User_function*)objectMedoidDistanceMin,
			false /* commute */, &op ), "MPI_Op_create for CObjectMedoidDistance" );
	}
	return op;
}

void MPIAPI CObjectMedoidDistance::objectMedoidDistanceMin(
	CObjectMedoidDistance* in, CObjectMedoidDistance* inout,
	int* length, MPI_Datatype* /*type*/ )
{
	for( int i = 0; i < *length; i++ ) {
		inout[i].Min( in[i] );
	}
}

///////////////////////////////////////////////////////////////////////////////

// Initializing or Build  step
void DoBuildStep( PamType& pam, CObjectMedoidDistance& best,
	const size_t objectBegin, const size_t objectEnd )
{
	best.Distance = numeric_limits<DistanceType>::max();
	best.Object = objectBegin;
	for( size_t object = objectBegin; object < objectEnd; object++ ) {
		if( pam.IsMedoid( object ) ) {
			continue; // if object is medoid
		}

		const DistanceType distance = ( pam.State() == PamType::Initializing ) ?
			pam.FindObjectDistanceToAll( object ) : -pam.AddMedoidProfit( object );

		if( distance < best.Distance ) {
			best.Distance = distance;
			best.Object = object;
		}
	}
}

// Swap step
void DoSwapStep( PamType& pam, CObjectMedoidDistance& best,
	const size_t objectBegin, const size_t objectEnd )
{
	best.Distance = 0;
	best.Medoid = pam.Medoids().front();
	best.Object = objectBegin;

	for( size_t object = objectBegin; object < objectEnd; object++ ) {
		if( pam.IsMedoid( object ) ) {
			continue; // if object is medoid
		}
		for( const size_t medoid : pam.Medoids() ) {
			DistanceType distance = pam.SwapResult( medoid, object );
			if( distance < best.Distance ) {
				best.Distance = distance;
				best.Medoid = medoid;
				best.Object = object;
			}
		}
	}
}

void CalcBeginEndObjects(
	const size_t numberOfObjects, const size_t numberOfProcess, const size_t rank,
	size_t& beginObject, size_t& endObject )
{
	const size_t objectsPerProcess = numberOfObjects / numberOfProcess;
	const size_t additionalObjects = numberOfObjects % numberOfProcess;
	beginObject = objectsPerProcess * rank + min( rank, additionalObjects );
	endObject = beginObject + objectsPerProcess;
	if( rank < additionalObjects ) {
		endObject++;
	}
}

void PamThread( PamType& pam,
	vector<CObjectMedoidDistance>& bests, CBarrier& barrier,
	size_t threadIndex,
	const size_t objectBegin, const size_t objectEnd )
{
#ifdef _DEBUG
	static mutex coutMutex;
#endif


	// Building and Initializing
	for( size_t i = 0; i < pam.NumberOfClusters(); i++ ) {
#ifdef _DEBUG
		{
			unique_lock<mutex> lock{ coutMutex };
			cout << CMpiSupport::Rank() << "," << threadIndex << " ["
				<< objectBegin << ", " << objectEnd << ") "
				<< "Building..." << i << endl;
		}
#endif
		DoBuildStep( pam, bests[threadIndex], objectBegin, objectEnd );

		barrier.Sync();

		if( threadIndex == 0 ) {
			for( CObjectMedoidDistance& objectMedoidDistance : bests ) {
				bests.front().Min( objectMedoidDistance );
			}
			bests.front().AllReduce();
			pam.AddMedoid( bests.front().Object );
		}

		barrier.Sync();
	}

	// Swapping
	for( size_t iteration = 0; iteration < 1000; iteration++ ) {
#ifdef _DEBUG
		{
			unique_lock<mutex> lock{ coutMutex };
			cout << CMpiSupport::Rank() << "," << threadIndex <<
				": " << "Swapping..." << iteration << endl;
		}
#endif
		DoSwapStep( pam, bests[threadIndex], objectBegin, objectEnd );

		barrier.Sync();

		if( threadIndex == 0 ) {
			for( CObjectMedoidDistance& objectMedoidDistance : bests ) {
				bests.front().Min( objectMedoidDistance );
			}
			bests.front().AllReduce();
			if( bests.front().Distance < 0 ) {
				pam.Swap( bests.front().Medoid, bests.front().Object );
			}
		}

		barrier.Sync();

		if( !( bests.front().Distance < 0 ) ) {
			break;
		}

		barrier.Sync();
	}
}

void DoPam( const size_t numberOfClusters, const DissimilarityMatrixType& matrix,
	size_t numberOfThreads = 1 )
{
	typedef CPartitioningAroundMedois<DissimilarityMatrixType> PamType;
	PamType pam( matrix, numberOfClusters );

	vector<thread> threads;
	threads.reserve( numberOfThreads );
	vector<CObjectMedoidDistance> bests( numberOfThreads );
	CBarrier barrier( numberOfThreads );

	for( size_t threadIndex = 0; threadIndex < numberOfThreads; threadIndex++ ) {
		size_t objectBegin = 0;
		size_t objectEnd = 0;
		CalcBeginEndObjects( pam.NumberOfObjects(),
			CMpiSupport::NumberOfProccess() * numberOfThreads,
			CMpiSupport::Rank() * numberOfThreads + threadIndex,
			objectBegin, objectEnd );

		threads.emplace_back( PamThread,
			ref( pam ), ref( bests ), ref( barrier ),
			threadIndex, objectBegin, objectEnd );
	}

	for( thread& t : threads ) {
		t.join();
	}

#ifdef _DEBUG
	if( CMpiSupport::Rank() == 0 ) {
		cout << endl;
		unordered_map<size_t, size_t> medoidToClusterId;
		for( size_t object = 0; object < pam.NumberOfObjects(); object++ ) {
			auto pair = medoidToClusterId.insert(
				make_pair( pam.ObjectMedoids()[object], medoidToClusterId.size() ) );
			cout << object << "\t" << pair.first->second << endl;
		}
	}
#endif
}

DissimilarityMatrixType BuildDissimilarityMatrixType( istream& input )
{
	size_t unused;
	size_t numberOfVectors;
	input >> unused >> numberOfVectors;
	CVector vector;
	CDissimilarityMatrixBuilder<CVector> builder;
	for( size_t i = 0; input.good() && i < numberOfVectors; i++ ) {
		input >> unused >> vector.X >> vector.Y;
		builder.push_back( vector );
	}
	if( !input.fail() ) {
		return builder.Build();
	}
	throw exception( "bad vectors file format!" );
}

void DoMain( const int argc, const char* const argv[] )
{
	if( argc < 3 || argc > 4 ) {
		throw exception( "too few arguments!\n"
			"Usage: pam NUMBER_OF_CLUSTERS VECTORS_FILENAME [NUMBER_OF_THREADS]" );
	}

	double readDataTime = 0.0;
	double pamTime = 0.0;

	const size_t numberOfClusters = stoul( argv[1] );
	DissimilarityMatrixType matrix;
	{
		CMpiTimer timer( readDataTime );
		matrix = BuildDissimilarityMatrixType( ifstream( argv[2] ) );
	}
	const size_t numberOfThreads = ( argc == 4 ) ? stoul( argv[3] ) : 1;

	{
		CMpiTimer timer( pamTime );
		DoPam( numberOfClusters, matrix, numberOfThreads );
	}

	cout << CMpiSupport::Rank() << "\t" << readDataTime << "\t" << pamTime << endl;
}

int main( int argc, char** argv )
{
	try {
		CMpiSupport::Initialize( &argc, &argv );
		DoMain( argc, argv );
		CMpiSupport::Finalize();
	} catch( exception& e ) {
		cerr << "Error: " << e.what() << endl;
		CMpiSupport::Abort( 1 );
		return 1;
	} catch( ... ) {
		cerr << "Unknown error!" << endl;
		CMpiSupport::Abort( 2 );
		return 2;
	}

	return 0;
}
