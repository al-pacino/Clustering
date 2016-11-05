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
	static size_t Rank();
	static size_t NumberOfProccess();

private:
	static bool initialized;
	static size_t rank;
	static size_t numberOfProccess;

	static void checkInitialized();
};

///////////////////////////////////////////////////////////////////////////////

inline size_t CMpiInitializer::Rank()
{
	checkInitialized();
	return rank;
}

inline size_t CMpiInitializer::NumberOfProccess()
{
	checkInitialized();
	return numberOfProccess;
}

///////////////////////////////////////////////////////////////////////////////
