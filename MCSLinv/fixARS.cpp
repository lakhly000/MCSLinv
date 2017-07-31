#include "fixARS.h"

/* Written by Richelle Streater, June 2017. */

/* FixARS regroups each storage vector into central angles,
to match with typical experimental data. It also divides each element
by its solid angle and number of particles and flips the ARS vector
so that normal transmittance is at 180 degrees rather than 0 degrees. */

/******************************************************************************/

void fixARS( vector<vector<vector<double> > >& ars, int numParticles,
        unsigned int m, unsigned int n ) {
	const double TAU = 6.28318530717958647692;
    double storeElement;
    double angleSize = 180.0 / ars.at(0).at(0).size();

    /* For loop applies these changes to the ARS vector for each combo of mut, etaa. */
    for ( unsigned int y = 0; y < m; y++ ) {
        for ( unsigned int z = 0; z < n; z++ ) {

            /* Add consecutive elements and store the results in the first half of ARS */
            for ( unsigned int i = 1; i < ars.at(y).at(z).size()/2; i++ ) {
                ars.at(y).at(z).at(i) = ars.at(y).at(z).at(2*i - 1)
                + ars.at(y).at(z).at(2*i);
            }

            /* Delete all elements over 1/2 of original size */
            ars.at(y).at(z).resize( ars.at(y).at(z).size() / 2 );

            /* Adjust by solid angle and number of particles:

                Solid angle for a sphere between two polar angles:
                Omega = 2*pi*(cos(theta1)-cos(theta2)

            */
            for( unsigned int j = 0; j < ars.at(y).at(z).size()-1; j++ ){
                ars.at(y).at(z).at( j+1 ) /= ( numParticles * TAU * ( cos( ( angleSize * ( j + 0.5 ) )
                    * TAU/180 ) - cos( ( angleSize * ( j + 1.5 ) ) * TAU/180 ) ) );
            }

            /* Double the first element to estimate a central angle of zero (since we don't have a -1 element) */
            ars.at(y).at(z).at( 0 ) /= ( numParticles * TAU * ( 1 - cos( angleSize / 2 * TAU/180 ) ) );

            /* Flip the vector */
            for ( unsigned int k = 0; k < ars.at(y).at(z).size() / 2; k++ ) {
                storeElement = ars.at(y).at(z).at(k);
                ars.at(y).at(z).at(k) = ars.at(y).at(z).at(ars.at(y).at(z).size()-1-k);
                ars.at(y).at(z).at(ars.at(y).at(z).size()-1-k) = storeElement;
            }
        }
    }
}
