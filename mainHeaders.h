#include "addVec.h"
#include "boundary.h"
#include "detect.h"
#include "fileToVec.h"
#include "fixARS.h"
#include "initSPRNG.h"
#include "layer.h"
#include "particle.h"
#include "propagate.h"
#include "dataOut.h"
#include "scatter.h"
#include "scoreParam.h"
#include "setParameters.h"
#include "specularR.h"
#include "subFromMax.h"
#include "updateInterval.h"
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <string>
#include <math.h>
#include "omp.h"

using namespace std;
