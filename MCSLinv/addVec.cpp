#include "addVec.h"

/* Written by Anne-Michelle Lieberson, July 2017. */

/* AddVec adds two 3D vectors element-wise and saves the contents to the first vector. Main calls
addVec to add the ARS results from each processor together. */

/******************************************************************************/

void addVec( vector<vector<vector<double> > > &v1, const vector<vector<vector<double> > > &v2, unsigned int dim1, unsigned int dim2 ) {
    for ( unsigned int i = 0; i < dim1; i++ ) {
        for ( unsigned int j = 0; j < dim2; j++ ) {
            for ( unsigned int k = 0; k < v1.at(0).at(0).size(); k++ ) {
                v1.at(i).at(j).at(k)+=v2.at(i).at(j).at(k);
            }
        }
    }
}
