#include "grid_data.h"


GridData::GridData() :
   mDfltValue(0.0), mMax(0.0,0.0,0.0) 
{
}

GridData::GridData(const GridData& orig) :
   mDfltValue(orig.mDfltValue)
{
   mData = orig.mData;
   mMax = orig.mMax;
}

GridData::~GridData() 
{
}

std::vector<double>& GridData::data()
{
   return mData;
}

GridData& GridData::operator=(const GridData& orig)
{
   if (this == &orig)
   {
      return *this;
   }
   mDfltValue = orig.mDfltValue;
   mData = orig.mData;
   mMax = orig.mMax;
   return *this;
}

void GridData::initialize(double dfltValue)
{
   mDfltValue = dfltValue;
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridData::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

const double GridData::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value  

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

void GridData::getCell(const vec3& pt, int& i, int& j, int& k)  
{
   vec3 pos = worldToSelf(pt); 
   i = (int) (pos[0]/theCellSize);
   j = (int) (pos[1]/theCellSize);
   k = (int) (pos[2]/theCellSize);   
}

double GridData::interpolate(const vec3& pt)
{
	
	// TODO: Implement sharper cubic interpolation here.  
	vec3 pos = worldToSelf(pt); 

	int i = (int) (pos[0]/theCellSize);
	int j = (int) (pos[1]/theCellSize);       
	int k = (int) (pos[2]/theCellSize);

	double scale = 1.0/theCellSize; 

	double fractx = scale*(pos[0] - i*theCellSize);
	double fracty = scale*(pos[1] - j*theCellSize);
	double fractz = scale*(pos[2] - k*theCellSize);

	assert (fractx < 1.0 && fractx >= 0);
	assert (fracty < 1.0 && fracty >= 0);
	assert (fractz < 1.0 && fractz >= 0);

	//(i, k)
	double tmp1 = (*this)(i,j-1,k);  
	double tmp2 = (*this)(i,j,k);  
	double tmp3 = (*this)(i,j+1,k);
	double tmp4 = (*this)(i,j+2,k);

	double a1 = tmp2; 
	double b1 = (tmp3-tmp1)/2.0; 
	double c1 = (tmp4-tmp2)/2.0; 
	double d1 = tmp3-tmp2;                

	double ik = HERMITE (a1,b1,c1,d1,fracty); 

	//(i+1, k)
	double tmp5 = (*this)(i+1,j-1,k);
	double tmp6 = (*this)(i+1,j,k);
	double tmp7 = (*this)(i+1,j+1,k);
	double tmp8 = (*this)(i+1,j+2,k);

    double a2 = tmp6; 
	double b2 = (tmp7-tmp5)/2.0; 
	double c2 = (tmp8-tmp6)/2.0; 
	double d2 = tmp7-tmp6; 

	double ip1k = HERMITE (a2,b2,c2,d2,fracty); 

	// (i, k+1)
	double tmp9 = (*this)(i,j-1,k+1);
	double tmp10 = (*this)(i,j,k+1);
	double tmp11 = (*this)(i,j+1,k+1);
	double tmp12 = (*this)(i,j+2,k+1);

	double a3 = tmp10; 
	double b3 = (tmp11-tmp9)/2.0; 
	double c3 = (tmp12-tmp10)/2.0; 
	double d3 = tmp11-tmp10;

	double ikp1 = HERMITE (a3,b3,c3,d3,fracty); 

	// (i+1, k+1)    
	double tmp13 = (*this)(i+1,j-1,k+1);
	double tmp14 = (*this)(i+1,j,k+1);
	double tmp15 = (*this)(i+1,j+1,k+1);
	double tmp16 = (*this)(i+1,j+2,k+1);

	double a4 = tmp14; 
	double b4 = (tmp15-tmp13)/2.0; 
	double c4 = (tmp16-tmp14)/2.0; 
	double d4 = (tmp15-tmp14); 

	double ip1kp1 = HERMITE (a4,b4,c4,d4,fracty);

	// (i-1,k)
	double tmp17 = (*this)(i-1,j-1,k);  
	double tmp18 = (*this)(i-1,j,k);  
	double tmp19 = (*this)(i-1,j+1,k);
	double tmp20 = (*this)(i-1,j+2,k);

	double a5 = tmp18; 
	double b5 = (tmp19 - tmp17)/2.0; 
	double c5 = (tmp20 - tmp18)/2.0; 
	double d5 = (tmp19 - tmp18);

	double is1k = HERMITE (a5,b5,c5,d5,fracty);

	//(i+2,k)
	double tmp21 = (*this)(i+2,j-1,k);  
	double tmp22 = (*this)(i+2,j,k);  
	double tmp23 = (*this)(i+2,j+1,k);
	double tmp24 = (*this)(i+2,j+2,k);

	double a6 = tmp22; 
	double b6 = (tmp23 - tmp21)/2.0; 
	double c6 = (tmp24 - tmp22)/2.0; 
	double d6 = (tmp23 - tmp22);   

	double ip2k = HERMITE (a6,b6,c6,d6,fracty);

	//(i-1, k+1)
	double tmp25 = (*this)(i-1,j-1,k+1);  
	double tmp26 = (*this)(i-1,j,k+1);  
	double tmp27 = (*this)(i-1,j+1,k+1);
	double tmp28 = (*this)(i-1,j+2,k+1);

	double a7 = tmp26; 
	double b7 = (tmp27 - tmp25)/2.0; 
	double c7 = (tmp28 - tmp26)/2.0; 
	double d7 = (tmp27 - tmp26);        

	double is1kp1 = HERMITE (a7,b7,c7,d7,fracty);


	//(i+2, k+1)
	double tmp29 = (*this)(i+2,j-1,k+1);  
	double tmp30 = (*this)(i+2,j,k+1);  
	double tmp31 = (*this)(i+2,j+1,k+1);
	double tmp32 = (*this)(i+2,j+2,k+1);

	double a8 = tmp30; 
	double b8 = (tmp31 - tmp29)/2.0; 
	double c8 = (tmp32 - tmp30)/2.0; 
	double d8 = (tmp31 - tmp30);  

	double ip2kp1 = HERMITE (a8,b8,c8,d8,fracty);

    double aK1 = ik;
	double bK1 = (ip1k - is1k)/2.0; 
	double cK1 = (ip2k - ik)/2.0; 
	double dK1 = ip1k- ik;

	double tmpK1 = HERMITE (aK1,bK1,cK1,dK1,fractx);

	double aK2 = ikp1;
	double bK2 = (ip1kp1 - is1kp1)/2.0; 
	double cK2 = (ip2kp1 - ikp1)/2.0; 
	double dK2 = ip1kp1 - ikp1;

	double tmpK2 = HERMITE (aK2,bK2,cK2,dK2,fractx);

	//(i-1, k-1) 
	double tmp33 = (*this)(i-1,j-1,k-1);  
	double tmp34 = (*this)(i-1,j,k-1);  
	double tmp35 = (*this)(i-1,j+1,k-1);
	double tmp36 = (*this)(i-1,j+2,k-1);

	double a9 = tmp34; 
	double b9 = (tmp35 - tmp33)/2.0; 
	double c9 = (tmp36 - tmp34)/2.0; 
	double d9 = tmp35 - tmp34; 

	double is1ks1 = HERMITE (a9,b9,c9,d9,fracty);

	//(i, k-1)
	double tmp37 = (*this)(i,j-1,k-1);  
	double tmp38 = (*this)(i,j,k-1);  
	double tmp39 = (*this)(i,j+1,k-1);
	double tmp40 = (*this)(i,j+2,k-1);

	double a10 = tmp38; 
	double b10 = (tmp39 - tmp37)/2.0; 
	double c10 = (tmp40 - tmp38)/2.0; 
	double d10 = tmp39 - tmp38; 

	double iks1 = HERMITE (a10,b10,c10,d10,fracty);
	
	//(i+1, k-1)
	double tmp41 = (*this)(i+1,j-1,k-1);  
	double tmp42 = (*this)(i+1,j,k-1);  
	double tmp43 = (*this)(i+1,j+1,k-1);
	double tmp44 = (*this)(i+1,j+2,k-1);

	double a11 = tmp42; 
	double b11 = (tmp43 - tmp41)/2.0; 
	double c11 = (tmp44 - tmp42)/2.0; 
	double d11 = tmp43 - tmp42; 

	double ip1ks1 = HERMITE (a11,b11,c11,d11,fracty);
    
	//(i+2, k-1)
	double tmp45 = (*this)(i+2,j-1,k-1);  
	double tmp46 = (*this)(i+2,j,k-1);  
	double tmp47 = (*this)(i+2,j+1,k-1);
	double tmp48 = (*this)(i+2,j+2,k-1);

	double a12 = tmp46; 
	double b12 = (tmp47 - tmp45)/2.0; 
	double c12 = (tmp48 - tmp46)/2.0; 
	double d12 = tmp47 - tmp46; 
	double ip2ks1 = HERMITE (a12,b12,c12,d12,fracty);


	//(i-1, k+2)
	double tmp49 = (*this)(i-1,j-1,k+2);  
	double tmp50 = (*this)(i-1,j,k+2);  
	double tmp51 = (*this)(i-1,j+1,k+2);
	double tmp52 = (*this)(i-1,j+2,k+2);

	double a13 = tmp50; 
	double b13 = (tmp51 - tmp49)/2.0; 
	double c13 = (tmp52 - tmp50)/2.0; 
	double d13 = tmp51 - tmp50; 
	double is1kp2 = HERMITE (a13,b13,c13,d13,fracty);

	//(i, k+2)
	double tmp53 = (*this)(i,j-1,k+2);  
	double tmp54 = (*this)(i,j,k+2);  
	double tmp55 = (*this)(i,j+1,k+2);
	double tmp56 = (*this)(i,j+2,k+2);

	double a14 = tmp54; 
	double b14 = (tmp55 - tmp53)/2.0; 
	double c14 = (tmp56 - tmp54)/2.0; 
	double d14 = tmp55 - tmp54;

	double ikp2 = HERMITE (a14,b14,c14,d14,fracty);

    //(i+1, k+2)
	double tmp57 = (*this)(i+1,j-1,k+2);  
	double tmp58 = (*this)(i+1,j,k+2);  
	double tmp59 = (*this)(i+1,j+1,k+2);
	double tmp60 = (*this)(i+1,j+2,k+2);

	double a15 = tmp58; 
	double b15 = (tmp59 - tmp57)/2.0; 
	double c15 = (tmp60 - tmp58)/2.0; 
	double d15 = tmp59 - tmp58;

	double ip1kp2 = HERMITE (a15,b15,c15,d15,fracty);

    //(i+2, k+2) 
	double tmp61 = (*this)(i+2,j-1,k+2);  
	double tmp62 = (*this)(i+2,j,k+2);  
	double tmp63 = (*this)(i+2,j+1,k+2);
	double tmp64 = (*this)(i+2,j+2,k+2);

	double a16 = tmp62; 
	double b16 = (tmp63 - tmp61)/2.0; 
	double c16 = (tmp64 - tmp62)/2.0; 
	double d16 = tmp63 - tmp62;

	double ip2kp2 = HERMITE (a16,b16,c16,d16,fracty);

	double aK0 = iks1;
	double bK0 = (ip1ks1 - is1ks1)/2.0; 
	double cK0 = (ip2ks1 - iks1)/2.0; 
	double dK0 = ip1ks1- iks1;

	double tmpK0 = HERMITE (aK0,bK0,cK0,dK0,fractx);

	double aK3 = ikp2;
	double bK3 = (ip1kp2 - is1kp2)/2.0; 
	double cK3 = (ip2kp2 - ikp2)/2.0; 
	double dK3 = ip1kp2- ikp2;

	double tmpK3 = HERMITE (aK3,bK3,cK3,dK3,fractx);

	double tFinal1 = tmpK1; 
	double tFinal2 = (tmpK2 - tmpK0)/2.0; 
	double tFinal3 = (tmpK3 - tmpK1)/2.0; 
	double tFinal4 = tmpK2 - tmpK1;

	double tmp = HERMITE (tFinal1,tFinal2,tFinal3,tFinal4,fractz);

	return tmp; 

	// LINEAR INTERPOLATION: 

    /*
	vec3 pos = worldToSelf(pt);

	int i = (int) (pos[0]/theCellSize);
	int j = (int) (pos[1]/theCellSize);       
	int k = (int) (pos[2]/theCellSize);

	double scale = 1.0/theCellSize;  
	double fractx = scale*(pos[0] - i*theCellSize);
	double fracty = scale*(pos[1] - j*theCellSize);
	double fractz = scale*(pos[2] - k*theCellSize);

	assert (fractx < 1.0 && fractx >= 0);
	assert (fracty < 1.0 && fracty >= 0);
	assert (fractz < 1.0 && fractz >= 0);

	// Y @ low X, low Z:
	double tmp1 = (*this)(i,j,k);   
	double tmp2 = (*this)(i,j+1,k);
	// Y @ high X, low Z:
	double tmp3 = (*this)(i+1,j,k);
	double tmp4 = (*this)(i+1,j+1,k);

	// Y @ low X, high Z:
	double tmp5 = (*this)(i,j,k+1);
	double tmp6 = (*this)(i,j+1,k+1);
	// Y @ high X, high Z:      
	double tmp7 = (*this)(i+1,j,k+1);
	double tmp8 = (*this)(i+1,j+1,k+1);

	// Y @ low X, low Z
	double tmp12 = LERP(tmp1, tmp2, fracty);
	// Y @ high X, low Z
	double tmp34 = LERP(tmp3, tmp4, fracty);

	// Y @ low X, high Z                             
	double tmp56 = LERP(tmp5, tmp6, fracty);
	// Y @ high X, high Z
	double tmp78 = LERP(tmp7, tmp8, fracty);

	// X @ low Z
	double tmp1234 = LERP (tmp12, tmp34, fractx);
	// X @ high Z
	double tmp5678 = LERP (tmp56, tmp78, fractx);

	// Z
	double tmp = LERP(tmp1234, tmp5678, fractz);   
	return tmp;*/

}

vec3 GridData::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0] - theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1] - theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2] - theCellSize*0.5), mMax[2]);
   return out; 
}

GridDataX::GridDataX() : GridData()
{
}

GridDataX::~GridDataX()
{
	
}

void GridDataX::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*(theDim[0]+1);
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize((theDim[0]+1)*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataX::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (i < 0 || i > theDim[0]) return dflt;

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData[stack + row + col];
}

const double GridDataX::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (i < 0 || i > theDim[0]) return dflt;

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData[stack + row + col];
}

vec3 GridDataX::worldToSelf(const vec3& pt) const
{   
   vec3 out;
   out[0] = min(max(0.0, pt[0]), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataY::GridDataY() : GridData()
{
}

GridDataY::~GridDataY()
{
}

void GridDataY::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*(theDim[1]+1);
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*(theDim[1]+1)*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataY::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (j < 0 || j > theDim[1]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData[stack + row + col];
}

const double GridDataY::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (j < 0 || j > theDim[1]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData[stack + row + col];
}

vec3 GridDataY::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataZ::GridDataZ() : GridData()
{
}

GridDataZ::~GridDataZ()
{
}

void GridDataZ::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*(theDim[2]+1);
   mData.resize(theDim[0]*theDim[1]*(theDim[2]+1), false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataZ::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (k < 0 || k > theDim[2]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData[stack + row + col];
}

const double GridDataZ::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (k < 0 || k > theDim[2]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData[stack + row + col];   
}

vec3 GridDataZ::worldToSelf(const vec3& pt) const    
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]), mMax[2]);
   return out;
}
