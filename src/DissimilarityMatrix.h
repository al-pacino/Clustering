#pragma once

///////////////////////////////////////////////////////////////////////////////

template<typename DISTANCE_TYPE>
class CDissimilarityMatrix {
	CDissimilarityMatrix( const CDissimilarityMatrix& ) = delete;
	CDissimilarityMatrix& operator=( const CDissimilarityMatrix& ) = delete;

public:
	typedef DISTANCE_TYPE DistanceType;

	CDissimilarityMatrix() :
		size( 0 )
	{
	}

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

	void Load( istream& input )
	{
		bool good = false;
		distances.clear();
		if( input.good() && input >> size ) {
			distances.reserve( size * size );
			DistanceType distance;
			while( input.good() && input >> distance ) {
				distances.push_back( distance );
			}
			if( distances.size() == size * size ) {
				good = true;
			}
		}
		if( !good ) {
			size = 0;
			distances.clear();
		}
	}

	void Save( ostream& output ) const
	{
		output << size;
		for( const DistanceType distance : distances ) {
			output << " " << distance;
		}
	}

protected:
	size_t size;
	vector<DistanceType> distances;
};

///////////////////////////////////////////////////////////////////////////////

template<typename OBJECT_TYPE>
class CDissimilarityMatrixBuilder :
	public vector<OBJECT_TYPE>,
	private CDissimilarityMatrix<typename OBJECT_TYPE::DistanceType> {
	CDissimilarityMatrixBuilder( const CDissimilarityMatrixBuilder& ) = delete;
	CDissimilarityMatrixBuilder& operator=( const CDissimilarityMatrixBuilder& ) = delete;

public:
	typedef OBJECT_TYPE ObjectType;
	typedef CDissimilarityMatrix<DistanceType> DissimilarityMatrixType;

	CDissimilarityMatrixBuilder()
	{
	}

	explicit CDissimilarityMatrixBuilder( size_t numberOfObjects )
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
