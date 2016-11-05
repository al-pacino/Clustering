#include <cmath>
#include <cassert>
#include <map>
#include <limits>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <exception>
#include <stdexcept>

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

struct CObjectMedoidDistance {
	unsigned long int Object;
	unsigned long int Medoid;
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
		MPI_Datatype types[count] = { MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, MPI_FLOAT };
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

		for( size_t i = 0; i < pam.Medoids().size(); i++ ) {
			const size_t medoid = pam.Medoids()[i];
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

void RunPam( PamType& pam, const size_t objectBegin, const size_t objectEnd )
{
	CObjectMedoidDistance best;

	// Building and Initializing
	for( size_t i = 0; i < pam.NumberOfClusters(); i++ ) {
#ifdef _DEBUG
		cout << CMpiSupport::Rank() << " ["
			<< objectBegin << ", " << objectEnd << ") "
			<< "Building..." << i << endl;
#endif
		DoBuildStep( pam, best, objectBegin, objectEnd );
		best.AllReduce();

		pam.AddMedoid( best.Object );
	}

	// Swapping
	for( size_t iteration = 0; iteration < 1000; iteration++ ) {
#ifdef _DEBUG
		cout << CMpiSupport::Rank()
			<< " Swapping..." << iteration << endl;
#endif
		DoSwapStep( pam, best, objectBegin, objectEnd );
		best.AllReduce();

		if( best.Distance < 0 ) {
			pam.Swap( best.Medoid, best.Object );
		} else {
			break;
		}
	}
}

void DoPam( const size_t numberOfClusters, const DissimilarityMatrixType& matrix )
{
	typedef CPartitioningAroundMedois<DissimilarityMatrixType> PamType;
	PamType pam( matrix, numberOfClusters );

	size_t objectBegin = 0;
	size_t objectEnd = 0;
	CalcBeginEndObjects( pam.NumberOfObjects(),
		CMpiSupport::NumberOfProccess(), CMpiSupport::Rank(),
		objectBegin, objectEnd );
	RunPam( pam, objectBegin, objectEnd );

#ifdef _DEBUG
	if( CMpiSupport::Rank() == 0 ) {
		cout << endl;
		map<size_t, size_t> medoidToClusterId;
		for( size_t object = 0; object < pam.NumberOfObjects(); object++ ) {
			pair<map<size_t, size_t>::iterator, bool> pair = medoidToClusterId.insert(
				make_pair( pam.ObjectMedoids()[object], medoidToClusterId.size() ) );
			cout << object << "\t" << pair.first->second << endl;
		}
		cout << endl;
	}
#endif
}

void BuildDissimilarityMatrix( istream& input, DissimilarityMatrixType& matrix )
{
	size_t unused = 0;
	size_t numberOfVectors = 0;
	input >> unused >> numberOfVectors;
	vector<CVector> vectors;
	vectors.resize( numberOfVectors );
	for( size_t i = 0; input.good() && i < numberOfVectors; i++ ) {
		input >> unused >> vectors[i].X >> vectors[i].Y;
	}
	if( input.fail() ) {
		throw domain_error( "bad vectors file format!" );
	}
	CDissimilarityMatrixBuilder<CVector>::Build( matrix, vectors.begin(), vectors.end() );
}

void DoMain( const int argc, const char* const argv[] )
{
	if( argc != 3 ) {
		throw domain_error( "too few arguments!\n"
			"Usage: pam NUMBER_OF_CLUSTERS VECTORS_FILENAME" );
	}

	double readDataTime = 0.0;
	double pamTime = 0.0;

	const size_t numberOfClusters = static_cast<size_t>( atoi( argv[1] ) );

	DissimilarityMatrixType matrix;
	{
		CMpiTimer timer( readDataTime );
		ifstream input( argv[2] );
		BuildDissimilarityMatrix( input, matrix );
	}

	{
		CMpiTimer timer( pamTime );
		DoPam( numberOfClusters, matrix );
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
