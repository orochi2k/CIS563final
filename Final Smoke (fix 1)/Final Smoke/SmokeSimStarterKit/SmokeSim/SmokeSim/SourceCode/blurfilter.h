#ifndef BLURFILTER_H
#define BLURFILTER_H

#include "grid_data.h"

class BlurFilter
{
public:
	BlurFilter();

	void SetSigma(double sigma);
	void GaussianBlur(const  GridData &src,  GridData &dst);

private:
	int mCenter;
	std::vector<double> mKernel;
	GridData mTemp;
};   


#endif // BLURFILTER_H    
