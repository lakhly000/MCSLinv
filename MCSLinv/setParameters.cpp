#include "setParameters.h"

/* Written by Richelle Streater, May 2017. */

/* SetParameters sets user-defined values from an input file. It is called
by main before the forward and inverse for loops. SetParameters defines a
stringstream to control which lines to include and which to ignore.
SetParameters ignores everything after a tab or a #, and completely ignores
contents of empty lines and lines that start with a #. This allows commenting
in the input file. For now, the material is considered to be single-layer.
SetParameters returns TRUE if there is no file error and FALSE if there is,
so that the program will close if there is a problem. */

/* Variables:
    layerVec: The layer object for the entire medium.
    mut, etaa: The vectors that hold the mut and etaa values to test for likelihood
    nVal, gVal, zMax: Parameters of the medium
    etaaMin, etaaMax, etaaN, mutMin, mutMax, mutN: Range and number of etaa and mut to test
    numParticles, numTrials: number of particles to start with (forward) and number of search
        iterations (inverse)
    numProc: Number of processors to use (parallelization)
    radius: Radius of detector (in experiment)

*/

/******************************************************************************/

bool setParameters( vector<Layer>& layerVec, vector<double>& mut, vector<double>& etaa,
    unsigned int& numParticles, unsigned int& numTrials, unsigned int& numProc, double& radius ){

    ifstream paramFile( "dataIn/input.txt" );
    unsigned int etaaN, mutN;
    double nVal, gVal, zMax, etaaMin, etaaMax, mutMin, mutMax;
    string line;
    stringstream l;

/*********************  Read file to stringstream  ****************************/

    if ( paramFile.is_open() ) {
        while( getline( paramFile, line ) ){

            /* Ignore lines that start with # or that have no contents */
                /* Ignore all line contents beyond tab or # */
            if ( ( line[0] != '#' ) && ( line.size()>0 ) ) {
                for (int i=0; i<line.size(); i++) {
                    if ( ( line[i] == '\t' ) || ( line[i] == '#' ) ) {
                        line.resize(i);
                        break;
                    }
                }
                l << line;
                l << '\n';
            }
        }

/*****************  End of reading file to stringstream  **********************/

/******************* Read stringstream to parameters  *************************/

        /* Read in all values from stringstream in order */
        l >> nVal;
        l >> gVal;
        l >> zMax;
        l >> etaaMin;
        l >> etaaMax;
        l >> etaaN;
        l >> mutMin;
        l >> mutMax;
        l >> mutN;
        l >> numParticles;
        l >> numTrials;
        l >> numProc;
        l >> radius;


        layerVec.at(0) = Layer( nVal, 1, 1, gVal, zMax );
        layerVec.at(0).setLayerNum(0);
    }

/***************  End of reading stringstream to parameters  ******************/

    else {
        cerr << "Error: File did not open (in setParameters.cpp)." << endl;
        return false;
    }

    /* Create etaa and mut vectors using defined N and range. */
    for ( unsigned int i = 0; i < etaaN; i++ ) {
        etaa.push_back( etaaMin + i * ( etaaMax-etaaMin ) / ( etaaN-1 ) );
    }

    for ( unsigned int i = 0; i < mutN; i++ ) {
        mut.push_back( mutMin + i * ( mutMax - mutMin ) / ( mutN - 1 ) );
    }

    paramFile.close();

    /* numProc is not defined when the file does not read properly */
    if ( numProc > 5000 ) {
        cerr << "File reading error (in setParameters.cpp)." << endl;
        return false;
    }
    return true;
}
