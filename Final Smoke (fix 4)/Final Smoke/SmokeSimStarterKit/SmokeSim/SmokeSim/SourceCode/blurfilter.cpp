#include "blurfilter.h"
#include "mac_grid.h"
#include "constants.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

BlurFilter::BlurFilter() :  
	mCenter(0)
{
}  

void BlurFilter::SetSigma(double sigma)
{
	assert(sigma > 0.0);

    int kernelSize = (int) (40 * sigma) * 2 + 1;
	mCenter = kernelSize / 2;
	
	//this is probably something I did wrong here 
	mKernel.resize(kernelSize);

    double normFactor = 1.0 / (sqrt(2.0 * M_PI) * sigma);

	for(int i=0; i < kernelSize; ++i)
	{
        int x = i - mCenter;
		mKernel[i] = normFactor * exp(-x * x / (2.0 * sigma * sigma));  
    }
}

void BlurFilter::GaussianBlur(const GridData &src, GridData &dst)   
{
    // dst = src;     
	//assert(src.NumRows() == dst.NumRows());
	//assert(src.NumCols() == dst.NumCols());
	//mTemp.Resize(src.NumRows(), src.NumCols());        
	 
	//Smooth horizontally    
    for(int j = 0; j < theDim[MACGrid::X]+2; j++)   //1)think about the limit here, it could be wrong; 2)also don't know if z is the axis we want to ignore.
	{                                               //there is something wrong with i,j. They don't change at all!!
        for(int i = 0; i < theDim[MACGrid::Y]+2; i++)  
		{
			double newValue = 0.0; 
			for(int cc = -mCenter; cc <= mCenter; ++cc)
			{
				if((((int)i + cc) >= 0) && (((int)i + cc) < (theDim[MACGrid::Y] +2)))
				{
					 newValue += src (i+cc, j, 0)* mKernel[mCenter+cc]; 
					//newValue += src (i+cc, j, 0); 
					//newValue += mKernel[mCenter+cc];          
				}
			}
			mTemp(i, j, 0) = newValue;    
		}
	}


	 // Smooth vertically.      
	for(int i=0; i < theDim[MACGrid::Y] + 2; ++i)
	{
		for(int j=0; j < theDim[MACGrid::X] + 2; ++j)
		{
            double newValue = 0.0;

			for(int rr = -mCenter; rr <= mCenter; ++rr)
			{
				if((((int)j + rr)>= 0) && (((int) j + rr) < (theDim[MACGrid::X]+2)))
				{
					newValue += mTemp(i, j+rr, 0) * mKernel[mCenter + rr];
				}
			}
			dst(i, j, 0) = newValue + 1e-256;  
		}
    }
}