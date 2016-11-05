#pragma once

#include <mpi.h>

///////////////////////////////////////////////////////////////////////////////

// Checking execution result of mpiFunctionName.
void MpiCheck( const int mpiResult, const string& mpiFunctionName );

///////////////////////////////////////////////////////////////////////////////

class CMpiSupport {
private:
	CMpiSupport();

public:
	static void Initialize( int* argc, char*** argv );
	static void Finalize();
	static void Abort( int code );
	static bool Initialized() { return initialized; }
	static size_t Rank();
	static size_t NumberOfProccess();
	static void Barrier();

private:
	static bool initialized;
	static size_t rank;
	static size_t numberOfProccess;

	static void checkInitialized();
};

///////////////////////////////////////////////////////////////////////////////

inline size_t CMpiSupport::Rank()
{
	checkInitialized();
	return rank;
}

inline size_t CMpiSupport::NumberOfProccess()
{
	checkInitialized();
	return numberOfProccess;
}

///////////////////////////////////////////////////////////////////////////////

class CMpiTimer {
private:
	CMpiTimer( const CMpiTimer& );
	CMpiTimer& operator=( const CMpiTimer& );

public:
	CMpiTimer( double& executionTime ) :
		time( executionTime ),
		startTime( getTime() )
	{
	}
	~CMpiTimer()
	{
		const double finishTime = getTime();
		time = finishTime - startTime;
	}

private:
	double& time;
	const double startTime;

	static double getTime()
	{
		CMpiSupport::Barrier();
		return MPI_Wtime();
	}
};

///////////////////////////////////////////////////////////////////////////////
