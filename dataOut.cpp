#include "dataOut.h"

/* Written by Richelle Streater, June 2017. */

/* DataOut creates the output file for the program. It creates a
file with the mean mua and mus values, the time-to-solution, and the
likelihood paraboloid coefficients for Mathematica to use. */

/******************************************************************************/

void dataOut( const vector<double>& dataList, int timeCount ) {
    string fnFolder = "dataOut/MCSLoutput.csv";

    /* Set up file to save to. */
    ofstream saveData;
    saveData.open( fnFolder.c_str() );

    /* Check for opening error and save file */
    if ( saveData.is_open() ) {
        saveData << "Seconds elapsed: " << timeCount << " s" << endl;
        saveData << "mu_s = " << dataList.at(4)*(1-dataList.at(3)) << endl;
        saveData << "mu_a = " << dataList.at(4)*dataList.at(3) << endl << endl;
        saveData << "Parameters for Mathematica: " << endl;
        for ( unsigned int i = 0; i < dataList.size()-1; i++ ) {
            saveData << dataList.at(i);
            saveData << endl;
        }
        saveData << dataList.back();
    }
    else {
        cerr << "File did not open (from dataOut.cpp)." << endl;
    }
}
