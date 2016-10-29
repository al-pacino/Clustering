#include <string>
#include <exception>

using namespace std;

#include <MpiInitializer.h>

///////////////////////////////////////////////////////////////////////////////

void MpiCheck( const int mpiResult, const string& mpiFunctionName )
{
	if( mpiResult != MPI_SUCCESS ) {
		throw exception( ( "MPI function '" + mpiFunctionName + "' failed"
			"with code '" + to_string( mpiResult ) + "'." ).c_str() );
	}
}

///////////////////////////////////////////////////////////////////////////////

bool CMpiInitializer::initialized = false;
int CMpiInitializer::rank = 0;
int CMpiInitializer::numberOfProccess = 0;

void CMpiInitializer::Initialize( int* argc, char*** argv )
{
	if( Initialized() ) {
		throw logic_error( "MPI was already initialized!" );
	}
	MpiCheck( MPI_Init( argc, argv ), "MPI_Init" );
	MpiCheck( MPI_Comm_rank( MPI_COMM_WORLD, &rank ), "MPI_Comm_rank" );
	MpiCheck( MPI_Comm_size( MPI_COMM_WORLD, &numberOfProccess ), "MPI_Comm_size" );
	initialized = true;
}

void CMpiInitializer::Finalize()
{
	checkInitialized();
	MPI_Finalize();
}

void CMpiInitializer::Abort( int code )
{
	if( Initialized() ) {
		MPI_Abort( MPI_COMM_WORLD, code );
	}
}

inline void CMpiInitializer::checkInitialized()
{
	if( !Initialized() ) {
		throw logic_error( "MPI was not initialized yet!" );
	}
}

///////////////////////////////////////////////////////////////////////////////
