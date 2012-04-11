// Written by Peter Kutz.

#ifndef CONSTANTS_H
#define CONSTANTS_H



#define LERP(a,b,t) (1-t)*a + t*b
#define HERMITE(a,b,c,d,t) a+b*t+(3*d-2*b-c)*t*t +(-2*d+b+c)*t*t*t; 


// Don't try to modify the values of these here.
// Modify the values of these in constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern const int theDim[3];
extern const double theCellSize;



#endif