// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#include <vector>
#undef max
#undef min
#include <fstream>


#include "blurfilter.h"
#define D_PRESSUREIT true
#define D_TEMPIT true
#define D_CLUTHU true
// Globals:
MACGrid target;

// NOTE: x -> cols, z -> rows, y -> stacks  
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
   for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 0; j < theDim[MACGrid::Y]; j++) \
         for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
   for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
      for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
         for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
         for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define FOR_EACH_CELLEX \
   for(int k = 1; k < theDim[MACGrid::Z]; k++)  \
      for(int j = 1; j < theDim[MACGrid::Y]; j++) \
         for(int i = 1; i < theDim[MACGrid::X]; i++) 



bool MACGrid::isValidFace(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] + 1 || j >= theDim[MACGrid::Y] + 1 || k >= theDim[MACGrid::Z] + 1) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}
	//if(this->mSolids(i,j,k)>0)
	//{
	//	return false;
	//}

	return true;
}

MACGrid::MACGrid()
{
   initialize();  

   //memory allocation for blurfilter pointer   
   mFilter = new BlurFilter;

   m_dx = 1.0 / theCellSize;
   m_2dx = 2.0 * m_dx;
   m_dt = 0.04; 
   m_dtdx = m_dt / m_dx;
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;    
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT; 

   //memory allocation for blurfilter pointer    
   mFilter = new BlurFilter;

   //set up dx here 
	m_dx = 1.0 / theCellSize;
	m_2dx = 2.0 * m_dx;
	m_dt = 0.04; 
   //set up dt/dx here 
	m_dtdx = m_dt / m_dx;


}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;  
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid()
{
	//release the blurfilter pointer 
	delete mFilter;
	mFilter = NULL;
}

void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);

   setUpAMatrix();  
}

void MACGrid::initialize()
{
   reset();
}

void MACGrid::updateSources()
{
    // TODO: Set initial values for density, temperature, and velocity.         
	mV (4,1,0) = 1; 
	mD (4,0,0) = 10.0; 
    mT (4,0,0) = 10;  
	/*mU(1,2,0) = 4;
	mD(0,2,0) = 0.8;
	mT(0,2,0) = 9;*/

	//for(int i=0; i < theDim[MACGrid::Y]; ++i)
	//{
	//		for(int j=0; j < theDim[MACGrid::X]; ++j)
	//		{
	//			//seems like I should not do anything about starting(source) density here, but the target density should be initialized later on
	//			//mDensity(i+1, j+1, 1) = srcDensity(i, j, 1);     
	//			mTargetDensity(i+1, j+1, 1) = tarDensity(i, j, 1);
	//		}
	//}

	//think about this where do you need to place the scaleMass ?  probably not 
	//scaleMass();
}

void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target.
	target.mU = mU;
    target.mV = mV;
    target.mW = mW;
	////////////////////////////////////////////////////q////////////////////////
	//vec3 totalforce(0,-9.8,0);
	FOR_EACH_FACE
	{
		//vec3 v = getVelocity(vec3(getCenter(i,j,k)));
		//vec3 poslast = vec3(i-v[0]*dt,j-v[1]*dt,k-v[2]*dt);
		//vec3 world = vec3(floor(i-v[0]*dt),floor(j-v[1]*dt),floor(k-v[2]*dt));
		//double aU = (poslast[0] - floor(poslast[0]));
		vec3 pnow = getCenter(i,j,k);
		//pnow -= vec3(theCellSize /2 , theCellSize /2,theCellSize/2);
		vec3 pnowX =  pnow - vec3(theCellSize /2 , 0,0);
		vec3 pnowY =  pnow - vec3(0 , theCellSize /2,0);
		vec3 pnowZ =  pnow - vec3(0 , 0,theCellSize /2);
		vec3 vX = getVelocity(vec3(pnowX));
		vec3 vY = getVelocity(vec3(pnowY));
		vec3 vZ = getVelocity(vec3(pnowZ));
		vec3 plastX = pnowX - vX * dt;
		vec3 plastY = pnowY - vY * dt;
		vec3 plastZ = pnowZ - vZ * dt;
		//if(plast[0] < 0) 
		//	plast[0] = 0;
		//if(plast[1] < 0) 
		//	plast[1] = 0;
		//if(plast[2] < 0) 
		//	plast[2] = 0;


		/*if(plastX[0] < 0) 
			plastX[0] = 0;
		if(plastX[1] < 0) 
			plastX[1] = 0;
		if(plastX[2] < 0) 
			plastX[2] = 0;
		if(plastY[0] < 0) 
			plastY[0] = 0;
		if(plastY[1] < 0) 
			plastY[1] = 0;
		if(plastY[2] < 0) 
			plastY[2] = 0;
		if(plastZ[0] < 0) 
			plastZ[0] = 0;
		if(plastZ[1] < 0) 
			plastZ[1] = 0;
		if(plastZ[2] < 0) 
			plastZ[2] = 0;*/


		//if(!this->areyousoild(plastX))
		target.mU(i,j,k) =  getVelocity(plastX)[0];
		//if(!this->areyousoild(plastY))
		target.mV(i,j,k) =  getVelocity(plastY)[1];
		//if(!this->areyousoild(plastZ))
		target.mW(i,j,k) =  getVelocity(plastZ)[2];
	}
	////////////////////////////////////////////////////////////////////////////
    // Then save the result to our object.
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target.
	target.mT = mT;
	FOR_EACH_CELL
	{
		vec3 pnow = getCenter(i,j,k);
		vec3 v = getVelocity(vec3(pnow));
		vec3 plast = pnow - v * dt;
		/*if(plast[0] < 0) 
			plast[0] = 0;
		if(plast[1] < 0) 
			plast[1] = 0;
		if(plast[2] < 0) 
			plast[2] = 0;*/
		/*if(!this->areyousoild(plast))*/
		target.mT(i,j,k) =  this->getTemperature(plast);
	}
    // Then save the result to our object.
    mT = target.mT;
}

void MACGrid::advectDensity(double dt)
{
    // TODO: Calculate new densitities and store in target.
	target.mD = mD;
	FOR_EACH_CELL
	{
		vec3 pnow = getCenter(i,j,k);
		vec3 v = getVelocity(vec3(pnow));
		vec3 plast = pnow - v * dt;
		/*if(plast[0] < 0) 
			plast[0] = 0;
		if(plast[1] < 0) 
			plast[1] = 0;
		if(plast[2] < 0) 
			plast[2] = 0;*/
		/*if(!this->areyousoild(plast))*/
		target.mD(i,j,k) =  this->getDensity(plast);
	}
    // Then save the result to our object.
    mD = target.mD;
}

void MACGrid::computeBouyancy(double dt)
{
	// TODO: Calculate bouyancy and store in target.
	target.mV = mV;
	bool bypass = false;
	if(D_TEMPIT)
	{
	///////////////////////////////////////////////////
		FOR_EACH_CELL
		{
			//target.mV(i,j,k) += dt * 0.5 * (target.mT(i,j,k) - 0);
			if (bypass || (isValidCell(i,j,k) && isValidCell(i,j+1,k)))
			{
			 target.mV(i,j+1,k) += dt * 0.5 * (target.mT(i,j,k)  - 0); 
			 //target.mV(i,j,k) += dt * 0.5 * (target.mT(i,j,k)  - 0) / 2; 
			}
		}
	///////////////////////////////////////////////////
	}
   // Then save the result to our object.
   mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt)
{
	// TODO: Calculate vorticity confinement forces.
	// Apply the forces to the current velocity and store the result in target.
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	////////////////////////
	GridData the_CLUTHU; 
	GridData the_FX; 
	GridData the_FY;
	GridData the_FZ; 
	the_FX.initialize(0.0);
	the_FY.initialize(0.0);
	the_FZ.initialize(0.0);
	the_CLUTHU.initialize(0.0);
	bool bypass = false;
	if( D_CLUTHU)
	{
       FOR_EACH_CELLEX
	   {
		   vec3 www = vec3();
		   www[0] =	(mW(i,j+1,k)-mW(i,j-1,k) - mV(i,j,k+1)+mV(i,j,k-1))/(2 *theCellSize);
		   www[1] = 	(mU(i,j,k+1)-mU(i,j,k-1) - mW(i+1,j,k)+mW(i-1,j,k))/(2*theCellSize);
		   www[2] =    (mV(i+1,j,k+1)-mV(i-1,j,k) - mU(i,j+1,k)+mU(i,j-1,k))/(2*theCellSize);
		   the_CLUTHU(i,j,k) = www.Length();
	   }

	 //  FOR_EACH_CELL
	 //  {
		//vec3 www = vec3();
		//www[0] =	(mW(i,j+1,k)-mW(i,j-1,k) - mV(i,j,k+1)+mV(i,j,k-1))/(2 *theCellSize);
		//www[1] = 	(mU(i,j,k+1)-mU(i,j,k-1) - mW(i+1,j,k)+mW(i-1,j,k))/(2*theCellSize);
		//www[2] =    (mV(i+1,j,k+1)-mV(i-1,j,k) - mU(i,j+1,k)+mU(i,j-1,k))/(2*theCellSize);
		//the_CLUTHU(i,j,k) = www.Length();
		///* vec3 pnow = getCenter(i,j,k);
		// vec3 v = getVelocity(vec3(pnow));
		// vec3 cluthu = vec3((v[2] - v[1])/theCellSize,(v[0]-v[2])/theCellSize,(v[1]-v[0])/theCellSize);
		// the_CLUTHU(i,j,k) = cluthu.Length();
		// if(i == 0 || j== 0||k== 0)
		// {
		//	 continue;
		// }
		// vec3 N;
		// N[0] = the_CLUTHU(i,j,k) - the_CLUTHU(i-1,j,k);
		// N[1] = the_CLUTHU(i,j,k) - the_CLUTHU(i,j-1,k);
		// N[2] = the_CLUTHU(i,j,k) - the_CLUTHU(i,j,k-1);*/
	 //  }
	      FOR_EACH_FACE
		 {
		 if(isValidFace(i,j,k) && isValidFace(i,j+1,k)&&isValidFace(i+1,j,k)&&isValidFace(i,j,k+1)&&isValidFace(i,j-1,k)&&isValidFace(i-1,j,k)&&isValidFace(i,j,k-1))
		 {
		vec3 www = vec3();
		www[0] =	(mW(i,j+1,k)-mW(i,j-1,k) - mV(i,j,k+1)+mV(i,j,k-1))/(2 *theCellSize);
		www[1] = 	(mU(i,j,k+1)-mU(i,j,k-1) - mW(i+1,j,k)+mW(i-1,j,k))/(2*theCellSize);
		www[2] =    (mV(i+1,j,k+1)-mV(i-1,j,k) - mU(i,j+1,k)+mU(i,j-1,k))/(2*theCellSize);

		 vec3 dwww = vec3();
		 dwww[0] = (the_CLUTHU(i+1,j,k)- the_CLUTHU(i-1,j,k))/ (2 * theCellSize);
		 dwww[1] = (the_CLUTHU(i,j+1,k)- the_CLUTHU(i,j-1,k))/ (2 * theCellSize);
		 dwww[2] = (the_CLUTHU(i,j,k+1)- the_CLUTHU(i,j,k-1))/ (2 * theCellSize);
		 vec3 nwww = dwww / (dwww.Length() + 0.0000000000000000001);
		 vec3 force = 0.8* theCellSize * (nwww.Cross(www));
		 if(i >= 1 && i <= theDim[MACGrid::X] - 1 )
		 {
		 target.mU(i,j,k) += force[0] * dt;
		 //target.mU(i,j,k) += force[0] * dt;
		 }
		// the_FX(i,j,k) = force[0] * dt;
		 if(j >= 1  && j <= theDim[MACGrid::Y] - 1)
		 {
		 target.mV(i,j,k) += force[1] * dt;
		 // target.mV(i,j,k) += force[1] * dt;
		 }
		// the_FY(i,j,k) = force[1] * dt;
		 if(k >= 1  && k <= theDim[MACGrid::Z] - 1)
		 {
		 target.mW(i,j,k) += force[2] * dt;
		 // target.mW(i,j,k) += force[2] * dt;
		 }
		// the_FZ(i,j,k) = force[2] * dt;
		 }
		 
	   }
	}
	////////////////////////
	// Then save the result to our object.
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
}
void MACGrid::addExternalForces(double dt)
{
     computeBouyancy(dt);
     computeVorticityConfinement(dt);
}

void MACGrid::project(double dt)
{
	// TODO: Solve Ap = d for pressure.
	// 1. Construct d
	// 2. Construct A
	// 3. Solve for p
	// Subtract pressure from our velocity and save in target.
	target.mP = mP;
	target.mU = mU;
	target.mV = mV;
	target.mW = mW;
	if(D_PRESSUREIT)
	{
	/////////////////////////////////////////////////////////////
	GridData the_d; 
	the_d.initialize(0.0);
	FOR_EACH_CELL
	{
		double vx1 = mU(i+1,j,k);
		double vx2 = mU(i,j,k);
		double vy1 = mV(i,j+1,k);
		double vy2 = mV(i,j,k);
		double vz1 = mW(i,j,k+1);
		double vz2 = mW(i,j,k);
		/*if(this->mSolids(i,j,k) > 0)
		{
			vx2 = 0;
			vy2 = 0;
			vz2 = 0;
		}*/
		/*if(this->mSolids(i+1,j,k) > 0)
		{
			vx1 = 0;
		}
		if(this->mSolids(i,j+1,k) > 0)
		{
			vy1 = 0;
		}
		if(this->mSolids(i,j,k+1) > 0)
		{
			vz1 = 0;
		}*/
		the_d(i,j,k) = (vx1 -vx2)/theCellSize + (vy1 -vy2)/theCellSize + (vz1 -vz2)/theCellSize;
		//double rou = target.mD(i,j,k) / 50.0 + 0.001;
		the_d(i,j,k) *= - (theCellSize * theCellSize) * 1 /dt;
		//the_d(i,j,k) *= -10;
	}
	setUpAMatrix();
	bool hit = conjugateGradient(this->AMatrix,target.mP,the_d,256,0.5);
	bool bypass = false;
		FOR_EACH_CELL
		{  
			/*if(mD(i,j,k) == 0)
			{
				continue;
			}*/
			//double change_P = (target.mP(i,j,k) - mP(i,j,k))/ theCellSize ;
			//change_P = target.mP(i,j,k) / theCellSize;
			//PRINT_LINE(change_P);
			//double rou = target.mD(i,j,k) / 50.0 + 0.001;
			double rou = 1;
			if (bypass || (isValidCell(i,j,k) && isValidCell(i-1,j,k)))
			{
			target.mU(i,j,k) = target.mU(i,j,k) - dt * 1.0 / rou * (target.mP(i,j,k) - target.mP(i-1,j,k))/ theCellSize;
			}
			if (bypass || (isValidCell(i,j,k) && isValidCell(i,j-1,k)))
			{
			target.mV(i,j,k) = target.mV(i,j,k) - dt * 1.0 / rou * (target.mP(i,j,k) - target.mP(i,j-1,k))/ theCellSize;
			}
			if (bypass || (isValidCell(i,j,k) && isValidCell(i,j,k-1))) 
			{
			target.mW(i,j,k) = target.mW(i,j,k) - dt * 1.0 / rou * (target.mP(i,j,k) - target.mP(i,j,k-1))/ theCellSize;
			}
		}
	/////////////////////////////////////////////////////////////
	}
	// Then save the result to our object
	mP = target.mP;
	mU = target.mU;
	mV = target.mV;
	mW = target.mW;
	// IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;       
   return vec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	} 

	return true;    
}

void MACGrid::setUpAMatrix() {

	FOR_EACH_CELL { 

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {   
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}

/////////////////////////////////////////////////////////////////////   

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}   

	GridData r = d; // Residual vector.

	GridData z; z.initialize();
	// TODO: Apply a preconditioner here.
	// For now, just bypass the preconditioner:   
	
	GridData precon; 
	GridData q; 
	precon.initialize(0);
	q.initialize(0); 

	double delta = 0.97;
    for(int i = 0; i < theDim[MACGrid::X]; i++) 
    {
	   for(int j = 0; j < theDim[MACGrid::Y]; j++)
	   {
		  for(int k = 0; k < theDim[MACGrid::Z]; k++)
			{
					double ep = A.diag (i,j,k) - (A.plusI (i-1,j,k)*precon (i-1,j,k))*(A.plusI (i-1,j,k)*precon (i-1,j,k))
												- (A.plusJ (i,j-1,k)*precon (i,j-1,k))*(A.plusJ (i,j-1,k)*precon (i,j-1,k))
												- (A.plusK (i,j,k-1)*precon (i,j,k-1))*(A.plusK (i,j,k-1)*precon (i,j,k-1))
												- delta*(A.plusI(i-1,j,k)*(A.plusJ(i-1,j,k) + A.plusK(i-1,j,k))*precon(i-1,j,k)*precon(i-1,j,k)
												+A.plusJ(i,j-1,k)*(A.plusI(i,j-1,k) + A.plusK(i,j-1,k))*precon(i,j-1,k)*precon(i,j-1,k)
												+A.plusK(i,j,k-1)*(A.plusI(i,j,k-1) + A.plusJ(i,j,k-1))*precon(i,j,k-1)*precon(i,j,k-1));
					precon (i,j,k) = 1.0/sqrt(ep); 
		 }  
	  }        
   }



	for(int i = 0; i < theDim[MACGrid::X]; i++) 
    {
	   for(int j = 0; j < theDim[MACGrid::Y]; j++)
	   {
		  for(int k = 0; k < theDim[MACGrid::Z]; k++)
	      {
				 double t1 = r(i,j,k) - A.plusI(i-1,j,k)*precon(i-1,j,k)*q(i-1,j,k)
									 - A.plusJ(i,j-1,k)*precon(i,j-1,k)*q(i,j-1,k)
									 - A.plusK(i,j,k-1)*precon(i,j,k-1)*q(i,j,k-1); 
				 q(i,j,k) = t1*precon(i,j,k); 
		 }
	  }        

   }

	for(int i = theDim[MACGrid::X]-1; i >= 0; i--) 
	{
		for(int j = theDim[MACGrid::Y]-1; j >= 0; j--) 
		{
	      for(int k = theDim[MACGrid::Z]-1; k >= 0; k--)  
		  {
			 
			
		    	 double t2 = q(i,j,k) - A.plusI(i,j,k)*precon(i,j,k)*z(i+1,j,k)
									 - A.plusJ(i,j,k)*precon(i,j,k)*z(i,j+1,k)
									 - A.plusK(i,j,k)*precon(i,j,k)*z(i,j,k+1); 
				 z(i,j,k) = t2*precon(i,j,k); 
			 
		 }  
	  }        
   }
 
	// z = r;
	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma;

		apply(A, s, z);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);

		if (maxMagnitude(r) <= tolerance) {

			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations."); 
			return true;
		}

		// TODO: Apply a preconditioner here.  
		// For now, just bypass the preconditioner:

	for(int i = 0; i < theDim[MACGrid::X]; i++) 
    {
	   for(int j = 0; j < theDim[MACGrid::Y]; j++)
	   {
		  for(int k = 0; k < theDim[MACGrid::Z]; k++)
	      {
			
				 double t1 = r(i,j,k) - A.plusI(i-1,j,k)*precon(i-1,j,k)*q(i-1,j,k)
									 - A.plusJ(i,j-1,k)*precon(i,j-1,k)*q(i,j-1,k)
									 - A.plusK(i,j,k-1)*precon(i,j,k-1)*q(i,j,k-1); 
				 q(i,j,k) = t1*precon(i,j,k); 
		 }
	  }        

   }

	for(int i = theDim[MACGrid::X]-1; i >= 0; i--) 
	{
		for(int j = theDim[MACGrid::Y]-1; j >= 0; j--) 
		{
	      for(int k = theDim[MACGrid::Z]-1; k >= 0; k--)  
		  {
		    	 double t2 = q(i,j,k) - A.plusI(i,j,k)*precon(i,j,k)*z(i+1,j,k)
									 - A.plusJ(i,j,k)*precon(i,j,k)*z(i,j+1,k)
									 - A.plusK(i,j,k)*precon(i,j,k)*z(i,j,k+1); 
				 z(i,j,k) = t2*precon(i,j,k); 
		 }
	  }        
   }
  

	   
	//   z = r; 

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;    

		sigma = sigmaNew;
	}

	//PRINT_LINE( "PCG didn't converge!" );
	return false;

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}

double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;    
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}


/////////////////////////////////////////////////////////////////////

void MACGrid::saveSmoke(const char* fileName) {
	std::ofstream fileOut(fileName);
	if (fileOut.is_open()) {
		FOR_EACH_CELL {
			fileOut << mD(i,j,k) << std::endl;
		}
		fileOut.close();
	}
}


//set up sigma to get prepared for driving force      
void MACGrid::SetParamSigma(double sigma)
{ 
	mFilter->SetSigma(sigma);

	//since our intial density are all zero, these codes should be all completely useless at this point 
    mFilter->GaussianBlur(mTargetDensity, mSmoothTargetDensity);   

	for(int i=1; i < theDim[MACGrid::Y]+2; ++i)
	{
        for(int j=1; j < theDim[MACGrid::X]+1; ++j)
		{
			double GradSmoothTargetDensity = (mSmoothTargetDensity(i, j, 1) -
												mSmoothTargetDensity(i-1, j, 1)) / m_dx;

			double AvgSmoothTargetDensity = (mSmoothTargetDensity(i, j, 1) +
												mSmoothTargetDensity(i-1, j, 1)) / 2.0;

            mNormalizedGradU (i, j, 1) = GradSmoothTargetDensity / AvgSmoothTargetDensity;
        }
    }

	for(int i=1; i < theDim[MACGrid::Y]+1; ++i)
	{
		for(int j=1; j < theDim[MACGrid::X]+2; ++j)
		{
			double GradSmoothTargetDensity = (mSmoothTargetDensity(i, j, 1) -
												mSmoothTargetDensity(i, j-1, 1)) / m_dx;

			double AvgSmoothTargetDensity = (mSmoothTargetDensity(i, j, 1) +
												mSmoothTargetDensity(i, j-1, 1)) / 2.0;

			mNormalizedGradV (i, j, 1) = GradSmoothTargetDensity / AvgSmoothTargetDensity;   
		}
    }    
}


void MACGrid::applyDrivingForce()
{
		mFilter->GaussianBlur(mD, mSmoothDensity);

        // Apply horizontal component of driving force.
		for(int i=1; i < theDim[MACGrid::Y] - 1; ++i)
		{
			for(int j=1; j < theDim[MACGrid::X] - 1; ++j)
			{
				double AvgSmoothDensity = (mSmoothDensity(i, j, 1) + mSmoothDensity(i-1, j, 1)) / 2.0;
				double Fu = AvgSmoothDensity * mNormalizedGradU(i, j, 1);

                mU(i, j, 1) += m_dt * m_vf * Fu;
			}
		}   

        // Apply vertical component of driving force.   
		for(int i=1; i < theDim[MACGrid::Y]- 1; ++i)
		{
			for(int j=1; j < theDim[MACGrid::X] - 1; ++j)
			{
				double AvgSmoothDensity = (mSmoothDensity(i, j, 1) + mSmoothDensity(i, j-1, 1)) / 2.0;
				double Fv = AvgSmoothDensity * mNormalizedGradV(i, j, 1);

                mV(i, j, 1) += m_dt * m_vf * Fv;   
            }
        }


		//need to set boundaries here. We are supposed to set all the boundaries to be 0 here 

		//setBoundaries(BND_H_VECTOR, mVelU);
        //setBoundaries(BND_V_VECTOR, mVelV);
}

void MACGrid::attenuateMomentum()
{
    // Attenuate horizontal component of momentum.
	for(int i=1; i < theDim[MACGrid::Y] - 1; ++i)
	{
		for(int j=1; j < theDim[MACGrid::X] - 1; ++j)
		{
            //mVelU(i, j) *= 1.0 - m_dt * m_vd;    n
                mU(i, j, 1) *= (1.0 - m_dt * m_vd);
		}
	}

	// Attenuate vertical component of momentum.    
	for(int i=1; i < theDim[MACGrid::Y] - 1; ++i)
	{
		for(int j=1; j < theDim[MACGrid::X] - 1; ++j)
		{
            //mVelV(i, j) *= 1.0 - m_dt * m_vd;
            mV(i, j, 1) *= (1.0 - m_dt * m_vd);
		}
	}


	//need to set boundaries here. We are supposed to set all the boundaries to be 0 here  


	//setBoundaries(BND_H_VECTOR, mVelU);
    //setBoundaries(BND_V_VECTOR, mVelV);
}



/////////////////////////////////////////////////////////////////////      

void MACGrid::draw(const Camera& c)
{   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities()
{
   // Draw line at each center
   //glColor4f(0.0, 1.0, 0.0, 1.0); // Use this if you want the lines to be a single color.
   glBegin(GL_LINES);
      FOR_EACH_CELL
      {
         vec3 pos = getCenter(i,j,k);
         vec3 vel = getVelocity(pos);
         if (vel.Length() > 0.0001)
         {
		   // Un-comment the line below if you want all of the velocity lines to be the same length.
           //vel.Normalize();  
           vel *= theCellSize/2.0; 
           vel += pos;
		   glColor4f(1.0, 1.0, 0.0, 1.0);
           glVertex3dv(pos.n);
		   glColor4f(0.0, 1.0, 0.0, 1.0);
           glVertex3dv(vel.n);
         }
      }
   glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k)
{

	// Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = mD(i, j, k); 
    return vec4(1.0, 0.0, 0.0, value);

}

vec4 MACGrid::getRenderColor(const vec3& pt)
{

	// TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
    double value = getDensity(pt); 
    return vec4(1.0, 0.0, 0.0, value);

}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front  
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front    
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}
