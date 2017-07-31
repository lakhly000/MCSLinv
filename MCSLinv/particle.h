#include "layer.h"
#include "weight.h"
#include <vector>

#ifdef SPRNGFIVE
#include "sprng_cpp.h"
#endif

using namespace std;

#pragma once

class Particle {
    public:
    #ifdef SPRNGFIVE
	Particle( double, vector<double>&, vector<double>&, Sprng* );
        Sprng* sprngptr;

    #else
    Particle( double, vector<double>&, vector<double>&, int* );
    int* sprngptr;
    #endif

	Particle();
	void reset( double );
	void updateWeightScatter();
	void updatePosition( double );
	vector<double> rVec;
	vector<double> dir;
    Layer lay;
    Weight weight;
};
