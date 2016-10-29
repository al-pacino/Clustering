#pragma once

#include <mpi.h>

///////////////////////////////////////////////////////////////////////////////

// Checking execution result of mpiFunctionName.
void MpiCheck( const int mpiResult, const string& mpiFunctionName );

///////////////////////////////////////////////////////////////////////////////

class CMpiInitializer {
	CMpiInitializer() = delete;

public:
	static void Initialize( int* argc, char*** argv );
	static void Finalize();
	static void Abort( int code );
	static bool Initialized() { return initialized; }
	static int Rank();
	static int NumberOfProccess();

private:
	static bool initialized;
	static int rank;
	static int numberOfProccess;

	static void checkInitialized();
};

///////////////////////////////////////////////////////////////////////////////

inline int CMpiInitializer::Rank()
{
	checkInitialized();
	return rank;
}

inline int CMpiInitializer::NumberOfProccess()
{
	checkInitialized();
	return numberOfProccess;
}

///////////////////////////////////////////////////////////////////////////////
