#include <cassert>
#include <limits>
#include <vector>
#include <iostream>
#include <exception>
#include <unordered_map>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

struct CPoint {
	typedef double NumbericType;
	typedef NumbericType DistanceType;

	NumbericType X;
	NumbericType Y;

	explicit CPoint( NumbericType x = 0, NumbericType y = 0 ) :
		X( x ), Y( y )
	{
	}

	DistanceType Distance( const CPoint& point )
	{
		const DistanceType dx = X - point.X;
		const DistanceType dy = Y - point.Y;
		return sqrt( dx * dx + dy * dy );
	}
};

///////////////////////////////////////////////////////////////////////////////

template<typename DISTANCE_TYPE>
class CDissimilarityMatrix {
	CDissimilarityMatrix( const CDissimilarityMatrix& ) = delete;
	CDissimilarityMatrix& operator=( const CDissimilarityMatrix& ) = delete;

public:
	typedef DISTANCE_TYPE DistanceType;

	CDissimilarityMatrix( CDissimilarityMatrix&& matrix )
	{
		*this = move( matrix );
	}

	CDissimilarityMatrix& operator=( CDissimilarityMatrix&& matrix )
	{
		size = matrix.size;
		matrix.size = 0;
		distances = move( matrix.distances );
		return *this;
	}

	size_t Size() const { return size; }

	DistanceType Distance( size_t i, size_t j ) const
	{
		assert( i < size && j < size );
		return distances[i * size + j];
	}

protected:
	size_t size;
	vector<DistanceType> distances;

	CDissimilarityMatrix()
	{
	}
};

///////////////////////////////////////////////////////////////////////////////

template<typename OBJECT_TYPE>
class CDissimilarityMatrixBuilder :
	public vector<OBJECT_TYPE>,
	private CDissimilarityMatrix<typename OBJECT_TYPE::DistanceType>
{
	CDissimilarityMatrixBuilder( const CDissimilarityMatrixBuilder& ) = delete;
	CDissimilarityMatrixBuilder& operator=( const CDissimilarityMatrixBuilder& ) = delete;

public:
	typedef OBJECT_TYPE ObjectType;
	typedef CDissimilarityMatrix<DistanceType> DissimilarityMatrixType;

	CDissimilarityMatrixBuilder()
	{
	}

	explicit CDissimilarityMatrixBuilder( size_t numberOfObject )
	{
		reserve( numberOfObjects );
	}

	CDissimilarityMatrix Build()
	{
		const size_t numberOfObjects = vector<ObjectType>::size();
		DissimilarityMatrixType::size = numberOfObjects;
		distances.reserve( numberOfObjects * numberOfObjects );
		for( size_t i = 0; i < numberOfObjects; i++ ) {
			for( size_t j = 0; j < numberOfObjects; j++ ) {
				if( i == j ) {
					distances.push_back( 0 );
				} else {
					distances.push_back( operator[]( i ).Distance( operator[]( j ) ) );
				}
			}
		}
		return move( static_cast<DissimilarityMatrixType&>( *this ) );
	}

	template<typename FORWARD_ITERATOR_TYPE>
	static DissimilarityMatrixType Build( FORWARD_ITERATOR_TYPE begin, FORWARD_ITERATOR_TYPE end )
	{
		CDissimilarityMatrixBuilder<OBJECT_TYPE> builder;
		while( begin != end ) {
			builder.push_back( *begin );
			++begin;
		}
		return builder.Build();
	}
};

///////////////////////////////////////////////////////////////////////////////

template<typename DISTNACE_TYPE>
class CPamClustering {
public:
	typedef DISTNACE_TYPE DistanceType;

	static void Pam( const CDissimilarityMatrix<DistanceType>& dissimilarityMatrix,
		size_t numberOfClusters, vector<size_t>& objectClusters )
	{
		objectClusters.clear();
		CPamClustering pam( dissimilarityMatrix, numberOfClusters );

		objectClusters.reserve( dissimilarityMatrix.Size() );
		unordered_map<size_t, size_t> medoidToIndex;
		for( const size_t objectMedoid : pam.objectMedoids ) {
			const size_t index = medoidToIndex.insert( make_pair( objectMedoid,
				medoidToIndex.size() ) ).first->second;
			objectClusters.push_back( index );
		}
	}

private:
	CPamClustering( const CPamClustering& ) = delete;
	CPamClustering& operator=( const CPamClustering& ) = delete;

	CPamClustering( const CDissimilarityMatrix<DistanceType>& dissimilarityMatrix,
			size_t numberOfClusters ) :
		matrix( dissimilarityMatrix )
	{
		if( numberOfClusters == 0 || numberOfClusters > dissimilarityMatrix.Size() ) {
			throw invalid_argument( "CPamClustering::CPamClustering invalid argument" );
		}

		if( numberOfClusters == 1 ) {
			// all objects are in same cluster
			objectMedoids.resize( dissimilarityMatrix.Size(), 0 );
		} else if( numberOfClusters == dissimilarityMatrix.Size() ) {
			// all objects are in own cluster
			objectMedoids.reserve( dissimilarityMatrix.Size() );
			for( size_t i = 0; i < dissimilarityMatrix.Size(); i++ ) {
				objectMedoids.push_back( i );
			}
		} else {
			build( numberOfClusters );
			swap();
		}
	}

	void build( size_t numberOfClusters );
	void swap();
	void calculateObjectMedoids();
	DistanceType swapMedoidAndObjectDistanceChange( size_t medoid,
		size_t j, size_t object ) const;

	bool isMedoid( size_t object ) const
	{
		return ( objectMedoids[object] == object );
	}

	DistanceType DistanceToMedoid( size_t object ) const
	{
		return matrix.Distance( object, objectMedoids[object] );
	}

	DistanceType DistanceToSecondMedoid( size_t object ) const
	{
		return matrix.Distance( object, objectSecondMedoids[object] );
	}

private:
	const CDissimilarityMatrix<DistanceType>& matrix;
	vector<size_t> medoids;
	vector<size_t> objectMedoids;
	vector<size_t> objectSecondMedoids;
};

///////////////////////////////////////////////////////////////////////////////

template<typename DISTNACE_TYPE>
void CPamClustering<DISTNACE_TYPE>::build( size_t numberOfClusters )
{
	medoids.reserve( numberOfClusters );

	DistanceType minDistance = numeric_limits<DistanceType>::max();
	size_t centralObject = 0;
	for( size_t i = 0; i < matrix.Size(); i++ ) {
		DistanceType distance = 0;
		for( size_t j = 0; j < matrix.Size(); j++ ) {
			distance += matrix.Distance( i, j );
		}
		if( distance < minDistance ) {
			minDistance = distance;
			centralObject = i;
		}
	}

	medoids.push_back( centralObject );
	objectMedoids.resize( matrix.Size(), centralObject );
	objectSecondMedoids.resize( matrix.Size() );

	for( size_t k = 1; k < numberOfClusters; k++ ) {
		DistanceType maxProfit = numeric_limits<DistanceType>::lowest();
		size_t bestObject = matrix.Size();
		for( size_t i = 0; i < matrix.Size(); i++ ) {
			if( isMedoid( i ) ) {
				continue; // if object is medoid
			}

			DistanceType profit = 0;
			for( size_t j = 0; j < matrix.Size(); j++ ) {
				if( i == j || isMedoid( j ) ) {
					continue; // if object is i or medoid
				}

				if( matrix.Distance( i, j ) < DistanceToMedoid( j ) ) {
					profit += DistanceToMedoid( j ) - matrix.Distance( i, j );
				}
			}
			if( maxProfit < profit ) {
				maxProfit = profit;
				bestObject = i;
			}
		}

		// add medoid
		medoids.push_back( bestObject );

		// calculate new objectMedoids
		for( size_t i = 0; i < matrix.Size(); i++ ) {
			if( isMedoid( i ) ) {
				continue; // if object is medoid
			}

			if( matrix.Distance( i, bestObject ) < DistanceToMedoid( i ) ) {
				objectMedoids[i] = bestObject;
			}
		}
	}
}

template<typename DISTNACE_TYPE>
void CPamClustering<DISTNACE_TYPE>::swap()
{
	while( true ) {
		calculateObjectMedoids();

		DistanceType minDistanceChange = 0;
		size_t bestMedoidIndex = medoids.size();
		size_t bestObject = matrix.Size();

		for( size_t medoidIndex = 0; medoidIndex < medoids.size(); medoidIndex++ ) {
			for( size_t object = 0; object < matrix.Size(); object++ ) {
				if( isMedoid( object ) ) {
					continue; // if object is medoid
				}

				DistanceType distanceChange = 0;
				for( size_t j = 0; j < matrix.Size(); j++ ) {
					if( j == object || isMedoid( j ) ) {
						continue; // if object is h or medoid
					}

					distanceChange += swapMedoidAndObjectDistanceChange(
						medoids[medoidIndex], j, object );
				}

				if( distanceChange < minDistanceChange ) {
					minDistanceChange = distanceChange;
					bestMedoidIndex = medoidIndex;
					bestObject = object;
				}
			}
		}

		if( minDistanceChange < 0 ) {
			medoids[bestMedoidIndex] = bestObject;
		} else {
			break;
		}
	}
}

template<typename DISTNACE_TYPE>
void CPamClustering<DISTNACE_TYPE>::calculateObjectMedoids()
{
	for( size_t i = 0; i < matrix.Size(); i++ ) {
		size_t objectMedoid = matrix.Size();
		DistanceType objectMedoidDistance = numeric_limits<DistanceType>::max();
		size_t objectSecondMedoid = matrix.Size();
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

		assert( objectMedoid < matrix.Size() && objectSecondMedoid < matrix.Size() );
		objectMedoids[i] = objectMedoid;
		objectSecondMedoids[i] = objectSecondMedoid;
	}
}

template<typename DISTNACE_TYPE>
DISTNACE_TYPE CPamClustering<DISTNACE_TYPE>::swapMedoidAndObjectDistanceChange(
	size_t medoid, size_t j, size_t object ) const
{
	if( objectMedoids[j] == medoid ) {
		// medoid is medoid of j-object
		if( DistanceToSecondMedoid( j ) > matrix.Distance( j, object ) ) {
			// object is new medoid of j-object
			return ( matrix.Distance( j, object ) - DistanceToMedoid( j ) );
		} else {
			// second j-object medoid is new medoid of j-object
			return ( DistanceToSecondMedoid( j ) - DistanceToMedoid( j ) );
		}
	} else {
		// medoid is NOT medoid of j-object
		if( DistanceToMedoid( j ) > matrix.Distance( j, object ) ) {
			// object is new medoid of j-object
			return ( matrix.Distance( j, object ) - DistanceToMedoid( j ) );
		} else {
			return 0;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
	try {

		vector<CPoint> points;
		points.emplace_back( 1.0, 1.0 );
		points.emplace_back( 2.0, 3.0 );
		points.emplace_back( 1.0, 2.0 );
		points.emplace_back( 2.0, 2.0 );
		points.emplace_back( 10.0, 4.0 );
		points.emplace_back( 11.0, 5.0 );
		points.emplace_back( 10.0, 6.0 );
		points.emplace_back( 12.0, 5.0 );
		points.emplace_back( 11.0, 6.0 );
		points.emplace_back( 5.0, 4.0 );
		points.emplace_back( 6.0, 3.0 );
		points.emplace_back( 6.0, 5.0 );
		points.emplace_back( 7.0, 4.0 );

		CDissimilarityMatrix<CPoint::DistanceType> matrix =
			CDissimilarityMatrixBuilder<CPoint>::Build( points.cbegin(), points.cend() );

		vector<size_t> pointClusters;
		CPamClustering<CPoint::DistanceType>::Pam( matrix, 3, pointClusters );

		for( size_t c : pointClusters ) {
			cout << c << endl;
		}
	} catch( exception& e ) {
		cerr << e.what() << endl;
		return 1;
	} catch( ... ) {
		cerr << "unknown error occured" << endl;
		return 1;
	}

	return 0;
}
