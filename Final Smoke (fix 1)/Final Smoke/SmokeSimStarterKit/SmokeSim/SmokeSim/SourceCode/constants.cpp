// Written by Peter Kutz.

#include "constants.h"

// TODO: Adjust the various settings below, and add more if you want.

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {4, 4, 4};    
#else
const int theDim[3] = {8, 8, 1};
#endif

const double theCellSize = 0.5;

