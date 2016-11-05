#pragma once

///////////////////////////////////////////////////////////////////////////////

template<typename DISTANCE_TYPE>
class CDissimilarityMatrix {
	template<typename OBJECT_TYPE> friend class CDissimilarityMatrixBuilder;

private:
	CDissimilarityMatrix( const CDissimilarityMatrix& matrix );
	CDissimilarityMatrix& operator=( const CDissimilarityMatrix& matrix );

public:
	typedef DISTANCE_TYPE DistanceType;

	CDissimilarityMatrix() :
		size( 0 )
	{
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
class CDissimilarityMatrixBuilder : public vector<OBJECT_TYPE> {
private:
	CDissimilarityMatrixBuilder( const CDissimilarityMatrixBuilder& );
	CDissimilarityMatrixBuilder& operator=( const CDissimilarityMatrixBuilder& );

public:
	typedef OBJECT_TYPE ObjectType;
	typedef typename OBJECT_TYPE::DistanceType DistanceType;
	typedef CDissimilarityMatrix<DistanceType> DissimilarityMatrixType;

	CDissimilarityMatrixBuilder( size_t numberOfObjects )
	{
		reserve( numberOfObjects );
	}

	void Build( DissimilarityMatrixType& matrix )
	{
		matrix.size = vector<ObjectType>::size();
		matrix.distances.clear();
		matrix.distances.reserve( matrix.size * matrix.size );
		for( size_t i = 0; i < matrix.size; i++ ) {
			for( size_t j = 0; j < matrix.size; j++ ) {
				if( i == j ) {
					matrix.distances.push_back( 0 );
				} else {
					matrix.distances.push_back( operator[]( i ).Distance( operator[]( j ) ) );
				}
			}
		}
	}

	template<typename FORWARD_ITERATOR_TYPE>
	static void Build( DissimilarityMatrixType& matrix, FORWARD_ITERATOR_TYPE begin, FORWARD_ITERATOR_TYPE end )
	{
		CDissimilarityMatrixBuilder<OBJECT_TYPE> builder( end - begin );
		while( begin != end ) {
			builder.push_back( *begin );
			++begin;
		}
		builder.Build( matrix );
	}
};

///////////////////////////////////////////////////////////////////////////////
