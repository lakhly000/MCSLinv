#include "mainHeaders.h"

/******************************  MCSLinv  *************************************/

/*
Program written by Richelle Streater, Anne-Michelle Lieberson, and Zachary Levine
July 2017
National Institute of Standards and Technology
Physical Measurement Lab
Sensor Science Division
Remote Sensing Group

Main first initializes all variables, seeds the RNG, and reads the input files. It then
zeros out the ars vectors (the solutions to the forward MC simulation) and sets the reference
importance sampling mua and mus values to the median of the input mut and etaa vectors. It
declares functions within the parallel for loop to eliminate race conditions. Next, it runs
the forward MC simulation on all particles, cycling between four states: scatter, propagate,
boundary, and detect, until the particle escapes or vanishes in the material. It updates the
search interval, doubles number of particles, and repeats again. */

/* Variables:
    numParticles: Number of particles to send through medium (forward MC simulation)
    numIter: Number of times to resize the search box (inverse)
    numProc: Number of processors to use (parallelization)
    seed, seedIn: seed for random number generator. If seedIn is 0, time(NULL) is used
    radius: Radius of detector
    layerVec: Vector holding parameters for all layers in the material
    mutVec, etaaVec: List of mut's and etaa's to search over (mut = mua+mus, etaa = mua/mut)
    expData: List of experimental angle resolved scattering results
    angleDiv: Number of divisions of ARS to measure
    T: Specular transmission, particles that initially make it into the medium
    ars, arsProc: Store the angle resolved scattering results (forward MC simulation)
    likGrid: The likelihood values for each point in the grid of mut and etaa
*/

int main() {

/**********************  Inputs and Initialization  ***************************/

    /* Read in experimental data, set parameters from file. */
    int seed, seedIn, time0 = time(NULL);
    unsigned int numParticles, numIter, numProc;
    double radius;
    vector<Layer> layerVec( 1 );
    Layer layAir;
    vector<double> mutVec, etaaVec, expData;
    fileToVec( expData, "dataIn/exp.txt" );

    /* Quit the program if there is an input error. */
    if ( !setParameters( layerVec, mutVec, etaaVec, numParticles,
        numIter, numProc, seedIn, radius ) ) {
        return 1;
    }

    if (seedIn == 0) {
       seed = time0;
    }

    else {
       seed = seedIn;
    }

    /* Seed the random number generator, depending on which version of
    SPRNG is present. */
#ifdef SPRNGFIVE
    int gtype = 2;
    /* Seed the RNG */
    Sprng *stream, **sprngptrarr;
    stream = SelectType(gtype);
    stream->init_sprng(0, 1, seed, SPRNG_DEFAULT);

    unsigned int nspawned = stream->spawn_sprng(numProc,&sprngptrarr);

    if(nspawned != numProc)
    {
        cout <<  "Error: " << nspawned << " actual out of " << numProc << " wanted streams spawned! " << endl;
        return 1;
    }

#else
    int **sprngptrarr;
    if ( !initSPRNG( numProc, &sprngptrarr, seed ) ) {
        return 1;
    }
#endif

    /* Initialize variables */

    /* Doubles, ints, and chars */
    unsigned int angleDiv = 2 * ( expData.size() );
    unsigned int mutSize = mutVec.size(), etaaSize = etaaVec.size();
    double T = 1 - specularR( layerVec.at(0) );

    /* Initialize vectors that need earlier information */
    omp_set_num_threads( numProc );
    vector<vector<vector<double> > > arsInitial( mutSize, vector<vector<double> >
        ( etaaSize, vector<double>( angleDiv, 0 ) ) ), ars;
    vector<vector<vector<vector<double> > > > arsProcInitial( numProc, vector<vector<vector<double> > >
        ( mutSize, vector<vector<double> >( etaaSize, vector<double>( angleDiv, 0 ) ) ) ), arsProc;
    vector<double> paramOut(5);
    vector<vector<double> > likGrid( mutSize, vector<double>( etaaSize, 0 ) );

    Layer *layPtr;
    layPtr = &layerVec.at(0);

/*********************  End of inputs and initialization  *********************/

/****************************  Inverse algorithm  *****************************/

    /* Control loop for inverse: resets bounding box and doubles particles each time. */
    for ( unsigned int a = 0; a < numIter; a++ ) {
        ars = arsInitial;
        arsProc = arsProcInitial;
        layerVec.at(0).setMua( mutVec.at( mutVec.size()/2 ) * etaaVec.at( etaaVec.size()/2 ) );
        layerVec.at(0).setMus( mutVec.at( mutVec.size()/2 ) - layerVec.at(0).getMua() );

/*********************  Forward Monte Carlo Simulation  ***********************/

        #pragma omp parallel for
        for ( unsigned int n = 0; n < numProc; n++ ) {

            /* Safe initializations within the parallel loop to eliminate race conditions */
            Particle par( T, mutVec, etaaVec, sprngptrarr[n] );
            int state;
            int propagate( Particle& );
            int detect( Particle&, double, unsigned int, vector<vector<vector<double> > >
                &, unsigned int, unsigned int );
            int scatter( Particle& );
            int boundary( Particle&, Layer&, vector<Layer>& );

            /* Send particle through the material */
            for ( unsigned int i = 0; i < numParticles/numProc; i++ ) {

                /* Reset particle weight/position and put the particle in the first layer */
                par.reset( T );
                par.lay = *layPtr;
                state = 2;

                /* Loop through states 1-4 until state is zero (AKA particle has escaped) */
                while ( state ) {
                    switch( state ) {
                    case 1:
                        state = scatter( par );
                        break;

                    case 2:
                        state = propagate( par );
                        break;

                    case 3:
                        state = boundary( par, layAir, layerVec );
                        break;

                    case 4:
                        state = detect( par, radius, angleDiv, arsProc.at(n), mutSize, etaaSize );
                        break;
                }
                }
            }
        }

        /* Add together parallel solutions to attain total ARS */
        for ( unsigned int p=0; p < numProc; p++ ) {
            addVec( ars, arsProc.at(p), mutSize, etaaSize );
        }

/******************  End of forward Monte Carlo simulation  *******************/

        /* Match ars format to experimental data, which shifts by half of an angle division. */
        fixARS( ars, numParticles, mutSize, etaaSize );

        /* Evaluate log-likelihood for each mut and etaa combination */
        if ( !scoreParam( ars, likGrid, expData, mutSize, etaaSize ) ) {
           return 1;
        }
        subFromMax( likGrid, mutSize, etaaSize );

        /* Resize the search region. Quit program if updateInterval has an error. */
        if ( !updateInterval( likGrid, paramOut, mutVec, etaaVec, a==(numIter-1) ) ) {
            return 1;
        }
        numParticles *= 2;
    }

/***********************  End of inverse algorithm ****************************/

    sprngptrarr = NULL;
    layPtr = NULL;

    dataOut( paramOut, time(NULL)-time0 );
    return 0;
}
