#include <cassert>
#include <limits>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <exception>
#include <algorithm>
#include <unordered_map>

using namespace std;

#include <MpiInitializer.h>
#include <Vector2d.h>
#include <DissimilarityMatrix.h>
#include <PartitioningAroundMedoids.h>

typedef float DistanceType;
typedef CVector2d<DistanceType> CVector;
typedef CDissimilarityMatrix<DistanceType> DissimilarityMatrixType;

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

	void AllReduce();

private:
	static MPI_Datatype datatype();
	static MPI_Op op();
	static void MPIAPI objectMedoidDistanceMin(
		CObjectMedoidDistance* in, CObjectMedoidDistance* inout,
		int* length, MPI_Datatype* /*type*/ );
};

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
		if( in[i].Distance < inout[i].Distance ) {
			inout[i] = in[i];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

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

void DoPam( const size_t numberOfClusters, const DissimilarityMatrixType& matrix )
{
	typedef CPartitioningAroundMedois<DissimilarityMatrixType> PamType;
	PamType pam( matrix, numberOfClusters );

	size_t objectBegin = 0;
	size_t objectEnd = 0;
	CalcBeginEndObjects( pam.NumberOfObjects(),
		CMpiInitializer::NumberOfProccess(), CMpiInitializer::Rank(),
		objectBegin, objectEnd );

	CObjectMedoidDistance best;

	// Building and Initializing
	for( size_t i = 0; i < pam.NumberOfClusters(); i++ ) {
#ifdef _DEBUG
		cout << "Building..." << i << endl;
#endif
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
		best.AllReduce();
		pam.AddMedoid( best.Object );
	}

	// Swapping
	for( size_t iteration = 0; iteration < 1000; iteration++ ) {
#ifdef _DEBUG
		cout << "Swapping..." << iteration << endl;
#endif
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

		best.AllReduce();

		if( best.Distance < 0 ) {
			pam.Swap( best.Medoid, best.Object );
		} else {
			break;
		}
	}

#ifdef _DEBUG
	if( CMpiInitializer::Rank() == 0 ) {
		unordered_map<size_t, size_t> medoidToClusterId;
		for( size_t object = 0; object < pam.NumberOfObjects(); object++ ) {
			auto pair = medoidToClusterId.insert(
				make_pair( pam.ObjectMedoids()[object], medoidToClusterId.size() ) );
			cout << object << "\t" << pair.first->second << endl;
		}
	}
#endif
}

void DoMain( const int argc, const char* const argv[] )
{
	if( argc != 3 ) {
		throw exception( "too few arguments!\n"
			"Usage: pam NUMBER_OF_CLUSTERS DISSIMILARITY_MATRIX_FILENAME" );
	}

	const size_t numberOfClusters = stoul( argv[1] );
	DissimilarityMatrixType matrix;
	matrix.Load( ifstream( argv[2] ) );

	DoPam( numberOfClusters, matrix );
}

int main( int argc, char** argv )
{
	try {
		CMpiInitializer::Initialize( &argc, &argv );
		DoMain( argc, argv );
		CMpiInitializer::Finalize();
	} catch( exception& e ) {
		cerr << "Error: " << e.what() << endl;
		CMpiInitializer::Abort( 1 );
		return 1;
	} catch( ... ) {
		cerr << "Unknown error!" << endl;
		CMpiInitializer::Abort( 2 );
		return 2;
	}

	return 0;
}
