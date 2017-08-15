#include "particle.h"

/* Particle is an object that represents a photon scattering in an optical medium. */

/* Members:
    rVec: (x,y,z) position coordinates
    dir: (kx,ky,kz) direction vector
    weight: an object carrying all weighting data for the particle (importance sampling)
    sprngptr: a pointer that eliminates race conditions in the parallel for loop
    layer: the layer that the particle is currently in
*/

/******************************************************************************/

#ifdef SPRNGFIVE
Particle::Particle ( double T, vector<double>& mutV, vector<double>& etaaV, Sprng* sprngptrin ) {

#else
Particle::Particle ( double T, vector<double>& mutV, vector<double>& etaaV, int* sprngptrin ) {
#endif

    /* Make a 3X1 position vector, filled with zeros */
    vector<double> positionVector( 3, 0 );

    /* Make 3X1 direction vector, point it down. */
    vector<double> dirVector( 3, 0 );
    dirVector.at( 2 ) = 1;

    rVec = positionVector;
    dir = dirVector;
    weight = Weight( mutV, etaaV );
    lay = Layer();
    sprngptr = sprngptrin;
}

Particle::Particle () {
    /* Make a 3X1 position vector, filled with zeros */
    vector<double> positionVector( 3, 0 );

    /* Make 3X1 direction vector, point it down. */
    vector<double> dirVector( 3, 0 );
    dirVector.at( 2 ) = 1;

    rVec = positionVector;
    dir = dirVector;
    lay = Layer();
    weight = Weight();
}

/* Reset function: creates new particle for beginning of main loop. */
void Particle::reset( double T ) {
    rVec = {0,0,0};
    dir = {0,0,1};
    weight.reset( T );
}

/* Update weight function- updates scalar w with MCML equation and updates
adjusted weight vector and weight vector for importance sampling. */
void Particle::updateWeightScatter() {
    weight.wScale *= lay.getMusDivMut();
    weight.updateWeightEtaa( lay.getDivVar() );
}

/* Update position function- Uses an input distance and particle direction. */
void Particle::updatePosition( double d ) {
    for ( unsigned int j = 0; j < rVec.size(); j++ ) {
        rVec.at(j) += d * dir.at(j);
    }
}
