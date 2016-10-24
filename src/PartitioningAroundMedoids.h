#pragma once

///////////////////////////////////////////////////////////////////////////////

template<typename DISSIMILARITY_MATRIX_TYPE>
class CPartitioningAroundMedois {
	CPartitioningAroundMedois( const CPartitioningAroundMedois& ) = delete;
	CPartitioningAroundMedois& operator=( const CPartitioningAroundMedois& ) = delete;

public:
	typedef DISSIMILARITY_MATRIX_TYPE DissimilarityMatrixType;
	typedef typename DissimilarityMatrixType::DistanceType DistanceType;

	enum StateType {
		Initializing,
		Building,
		Swapping
	};

	explicit CPartitioningAroundMedois( const DissimilarityMatrixType& dissimilarityMatrix,
			size_t _numberOfClusters ) :
		matrix( dissimilarityMatrix ),
		numberOfClusters( _numberOfClusters ),
		state( Initializing )
	{
		if( numberOfClusters < 2 || numberOfClusters > matrix.Size() ) {
			throw invalid_argument( "CPartitioningAroundMedois initializing failed" );
		}

		medoids.reserve( numberOfClusters );
		objectMedoids.resize( matrix.Size() );
		objectSecondMedoids.resize( matrix.Size() );
	}

	const DissimilarityMatrixType& DissimilarityMatrix() const { return matrix; }
	size_t NumberOfObjects() const { return matrix.Size(); }
	size_t NumberOfClusters() const { return numberOfClusters; }
	StateType State() const { return state; }
	const vector<size_t>& Medoids() const { return medoids; }
	const vector<size_t>& ObjectMedoids() const { return objectMedoids; }
	bool IsMedoid( size_t object ) const
	{
		return ( objectMedoids[object] == object );
	}
	// build operations
	DistanceType FindObjectDistanceToAll( size_t object ) const;
	void AddMedoid( size_t object );
	DistanceType AddMedoidProfit( size_t object ) const;
	// swap operatations
	void Swap( size_t medoid, size_t object );
	DistanceType SwapResult( size_t medoid, size_t object ) const;

private:
	const DissimilarityMatrixType& matrix;
	const size_t numberOfClusters;
	StateType state;
	vector<size_t> medoids;
	vector<size_t> objectMedoids;
	vector<size_t> objectSecondMedoids;

	DistanceType distanceToMedoid( size_t object ) const
	{
		return matrix.Distance( object, objectMedoids[object] );
	}
	DistanceType distanceToSecondMedoid( size_t object ) const
	{
		return matrix.Distance( object, objectSecondMedoids[object] );
	}
	void findObjectMedoids();
	DistanceType swapResult( size_t medoid, size_t j, size_t object ) const;
};

///////////////////////////////////////////////////////////////////////////////

template<typename DMT>
typename CPartitioningAroundMedois<DMT>::DistanceType
CPartitioningAroundMedois<DMT>::FindObjectDistanceToAll( size_t object ) const
{
	assert( object < NumberOfObjects() );

	DistanceType distance = 0;
	for( size_t anotherObject = 0; anotherObject < NumberOfObjects(); anotherObject++ ) {
		distance += matrix.Distance( object, anotherObject );
	}

	return distance;
}

template<typename DMT>
void CPartitioningAroundMedois<DMT>::AddMedoid( size_t medoid )
{
	assert( State() == Initializing || State() == Building );
	assert( medoid < NumberOfObjects() );
	assert( medoids.empty() == ( State() == Initializing ) );
	assert( medoids.size() < NumberOfClusters() );

	medoids.push_back( medoid );

	if( State() == Initializing ) {
		for( size_t& objectMedoid : objectMedoids ) {
			objectMedoid = 0;
		}

		state = Building;
	} else {
		// calculate new objectMedoids
		for( size_t object = 0; object < NumberOfObjects(); object++ ) {
			if( IsMedoid( object ) ) {
				continue; // if object is medoid
			}

			if( matrix.Distance( object, medoid ) < distanceToMedoid( object ) ) {
				objectMedoids[object] = medoid;
			}
		}

		if( medoids.size() == NumberOfClusters() ) {
			state = Swapping;

			findObjectMedoids();
		}
	}
}

template<typename DMT>
typename CPartitioningAroundMedois<DMT>::DistanceType
CPartitioningAroundMedois<DMT>::AddMedoidProfit( size_t object ) const
{
	assert( State() == Building );
	assert( object < NumberOfObjects() );
	assert( !IsMedoid( object ) );

	DistanceType profit = 0;
	for( size_t anotherObject = 0; anotherObject < NumberOfObjects(); anotherObject++ ) {
		if( object == anotherObject || IsMedoid( anotherObject ) ) {
			continue; // if anotherObject is object or medoid
		}

		if( matrix.Distance( object, anotherObject ) < distanceToMedoid( anotherObject ) ) {
			profit += distanceToMedoid( anotherObject ) - matrix.Distance( object, anotherObject );
		}
	}

	return profit;
}

template<typename DMT>
void CPartitioningAroundMedois<DMT>::Swap( size_t medoid, size_t object )
{
	assert( object < NumberOfObjects() );
	assert( State() == Swapping );

	auto mi = find( medoids.begin(), medoids.end(), medoid );
	assert( mi != medoids.end() );
	*mi = object;

	findObjectMedoids();
}

template<typename DMT>
typename CPartitioningAroundMedois<DMT>::DistanceType
CPartitioningAroundMedois<DMT>::SwapResult( size_t medoid, size_t object ) const
{
	assert( medoid < NumberOfObjects() );
	assert( object < NumberOfObjects() );
	assert( IsMedoid( medoid ) );
	assert( !IsMedoid( object ) );
	assert( State() == Swapping );

	DistanceType result = 0;
	for( size_t j = 0; j < NumberOfObjects(); j++ ) {
		if( j == object || IsMedoid( j ) ) {
			continue; // if j is object or medoid
		}
		result += swapResult( medoid, j, object );
	}
	return result;
}

template<typename DMT>
void CPartitioningAroundMedois<DMT>::findObjectMedoids()
{
	assert( State() != Initializing );

	for( size_t i = 0; i < NumberOfObjects(); i++ ) {
		size_t objectMedoid = NumberOfObjects();
		DistanceType objectMedoidDistance = numeric_limits<DistanceType>::max();
		size_t objectSecondMedoid = NumberOfObjects();
		DistanceType objectSecondMedoidDistance = numeric_limits<DistanceType>::max();

		for( const size_t medoid : medoids ) {
			const DistanceType distance = matrix.Distance( medoid, i );
			if( distance < objectMedoidDistance ) {
				objectSecondMedoid = objectMedoid;
				objectSecondMedoidDistance = objectMedoidDistance;
				objectMedoid = medoid;
				objectMedoidDistance = distance;
			} else if( distance < objectSecondMedoidDistance ) {
				objectSecondMedoid = medoid;
				objectSecondMedoidDistance = distance;
			}
		}

		assert( objectMedoid < NumberOfObjects() && objectSecondMedoid < NumberOfObjects() );
		objectMedoids[i] = objectMedoid;
		objectSecondMedoids[i] = objectSecondMedoid;
	}
}

template<typename DMT>
typename CPartitioningAroundMedois<DMT>::DistanceType
CPartitioningAroundMedois<DMT>::swapResult( size_t medoid, size_t j, size_t object ) const
{
	if( objectMedoids[j] == medoid ) {
		// medoid is medoid of j-object
		if( distanceToSecondMedoid( j ) > matrix.Distance( j, object ) ) {
			// object is new medoid of j-object
			return ( matrix.Distance( j, object ) - distanceToMedoid( j ) );
		} else {
			// second j-object medoid is new medoid of j-object
			return ( distanceToSecondMedoid( j ) - distanceToMedoid( j ) );
		}
	} else {
		// medoid is NOT medoid of j-object
		if( distanceToMedoid( j ) > matrix.Distance( j, object ) ) {
			// object is new medoid of j-object
			return ( matrix.Distance( j, object ) - distanceToMedoid( j ) );
		} else {
			return 0;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
