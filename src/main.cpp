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

int main( int argc, const char* argv[] )
{
	try {
		if( argc != 3 ) {
			throw invalid_argument( "too few arguments!\n"
				"Usage: pam NUMBER_OF_CLUSTERS DISSIMILARITY_MATRIX_FILENAME" );
		}

		typedef CVector2d<float> CVector;
		typedef CVector::DistanceType DistanceType;
		typedef CDissimilarityMatrix<DistanceType> DissimilarityMatrixType;

		const size_t numberOfClusters = stoul( argv[1] );
		DissimilarityMatrixType matrix;
		matrix.Load( ifstream( argv[2] ) );
		CPartitioningAroundMedois<DissimilarityMatrixType> pam( matrix, numberOfClusters );

		// Initializing
		cout << "Building...0" << endl;
		size_t centralObject = 0;
		CVector::DistanceType minDistance = pam.FindObjectDistanceToAll( 0 );
		for( size_t object = 1; object < pam.NumberOfObjects(); object++ ) {
			const CVector::DistanceType distance = pam.FindObjectDistanceToAll( object );
			if( distance < minDistance ) {
				minDistance = distance;
				centralObject = object;
			}
		}
		pam.AddMedoid( centralObject );

		// Building
		for( size_t i = 1; i < pam.NumberOfClusters(); i++ ) {
			cout << "Building..." << i << endl;

			size_t medoid = pam.NumberOfObjects();
			CVector::DistanceType maxDistance = 0;
			for( size_t object = 0; object < pam.NumberOfObjects(); object++ ) {
				if( pam.IsMedoid( object ) ) {
					continue; // if object is medoid
				}

				const CVector::DistanceType distance = pam.AddMedoidProfit( object );
				if( distance > maxDistance ) {
					maxDistance = distance;
					medoid = object;
				}
			}
			pam.AddMedoid( medoid );
		}

		cout << endl;

		// Swapping
		for( size_t iteration = 0; iteration < 1000; iteration++ ) {
			cout << "Swapping..." << iteration << endl;

			CVector::DistanceType bestResult = 0;
			size_t bestMedoid = pam.NumberOfObjects();
			size_t bestObject = pam.NumberOfObjects();

			for( const size_t medoid : pam.Medoids() ) {
				for( size_t object = 0; object < pam.NumberOfObjects(); object++ ) {
					if( pam.IsMedoid( object ) ) {
						continue; // if object is medoid
					}

					CVector::DistanceType result = pam.SwapResult( medoid, object );
					if( result < bestResult ) {
						bestResult = result;
						bestMedoid = medoid;
						bestObject = object;
					}
				}
			}

			if( bestResult < 0 ) {
				pam.Swap( bestMedoid, bestObject );
			} else {
				break;
			}
		}

		cout << endl;

		unordered_map<size_t, size_t> medoidToClusterId;
		for( size_t object = 0; object < pam.NumberOfObjects(); object++ ) {
			auto pair = medoidToClusterId.insert(
				make_pair( pam.ObjectMedoids()[object], medoidToClusterId.size() ) );
			cout << object << "\t" << pair.first->second << endl;
		}
	} catch( exception& e ) {
		cerr << "Error: " << e.what() << endl;
		return 1;
	}

	return 0;
}
