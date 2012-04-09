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

MACGrid::MACGrid()
{
   initialize();   
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;    
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT; 
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
	mU (4,7,0) = 0; 
	mV (4,7,0) = 10; 
	mW (4,7,0) = 0; 
	mD (4,0,0) = 10.0; 
    mT (4,0,0) = 10;  
}

void MACGrid::advectVelocity(double dt)
{
    // TODO: Calculate new velocities and store in target.          
   target.mV = mV;
   target.mU = mU;
   target.mW = mW;

   for(int k = 0; k < theDim[MACGrid::Z]; k++) 
   {
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; i++) 
		 {
	        //consider the velocity V first    
			double xAxis = (double)(theCellSize*(i+i+1+i+i+1)/4.0);
			double yAxis = (double)(theCellSize*j); 
		    double zAxis = (double)(theCellSize*(k+k+1+k+k+1)/4.0);
			vec3 pV  = vec3 (xAxis, yAxis, zAxis);   
			vec3 cVel = getVelocity (pV);        
			vec3 oldPV = pV -dt*cVel; 

			//check boundary for oldPV  
			if(oldPV[0] <= 0) 
			{
				oldPV[0] = 0; 
			}
			if(oldPV[1] <= 0) 
			{
				oldPV[1] = 0; 
			}
			if(oldPV[2] <= 0) 
			{
				oldPV[2] = 0; 
			}

			vec3 oldcVel = getVelocity (oldPV);  

			target.mU (i, j, k) = oldcVel[0]; 
			target.mV (i, j, k) = oldcVel[1]; 
			target.mW (i, j, k) = oldcVel[2]; 
		 }
	  }   
	}       

   //only update the x axis 
   for(int k = 0; k < theDim[MACGrid::Z]; k++) 
   {
      for(int j = 0; j <  theDim[MACGrid::Y]; j++) 
	  {
         for(int i = 0; i <  theDim[MACGrid::X]+1; i++) 
		 {
	        //consider the velocity V      I
		    double xAxis = (double) (theCellSize*i); 
			double yAxis = (double)(theCellSize*(j+j+1+j+j+1)/4.0); 
			double zAxis = (double)(theCellSize*(k+k+1+k+k+1)/4.0); 
			vec3 pU  = vec3 (xAxis, yAxis, zAxis);  
			vec3 cVel = getVelocity (pU);
			vec3 oldPU = pU -dt*cVel;     

			//check boundary for oldPV     
			if(oldPU[0] <= 0) 
			{
				oldPU[0] = 0; 
			}
			if(oldPU[1] <= 0)  
			{
				oldPU[1] = 0; 
			}
			if(oldPU[2] <= 0) 
			{
				oldPU[2] = 0; 
			}

			vec3 oldcVel = getVelocity (oldPU);

			target.mU (i, j, k) = oldcVel[0];
			target.mV (i, j, k) = oldcVel[1];
			target.mW (i, j, k) = oldcVel[2];
		 }
	  }  
	}

   //update in z axis   
   for(int k = 0; k < theDim[MACGrid::Z]+1; k++) 
   {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; i++) 
		 {
	        //consider the velocity V      
			double xAxis = (double)(theCellSize*(i+i+1+i+i+1)/4.0); 
			double yAxis = (double)(theCellSize*(j+j+1+j+j+1)/4.0); 
			double zAxis = (double)(theCellSize*k);   
			vec3 pW  = vec3 (xAxis, yAxis, zAxis);  
			vec3 cVel = getVelocity (pW); 
            vec3 oldPW = pW - dt*cVel; 

			//check boundary for oldPV    
			if(oldPW[0] <= 0) 
			{
				oldPW[0] = 0; 
			}
			if(oldPW[1] <= 0) 
			{
				oldPW[1] = 0; 
			}
			if(oldPW[2] <= 0) 
			{
				oldPW[2] = 0; 
			}
			vec3 oldcVel = getVelocity (oldPW);
			target.mU (i, j, k) = oldcVel[0];
			target.mV (i, j, k) = oldcVel[1];
			target.mW (i, j, k) = oldcVel[2];
		 }
	  }
	}
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectTemperature(double dt)
{
    // TODO: Calculate new temp and store in target.
   target.mT = mT;
   for(int k = 0; k < theDim[MACGrid::Z]; k++) 
   {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; i++) 
		 { 
			vec3 P = getCenter(i, j, k);
			vec3 temV = getVelocity(P);  
			vec3 oldP = P -dt*temV; 
			if(oldP[0] <= 0)    
			{
				oldP = vec3(0, oldP[1], oldP[2]); 
			}
			if(oldP[1] <= 0) 
			{
				oldP = vec3(oldP[0], 0, oldP[2]);   
			}
			if(oldP[2] <= 0) 
			{
				oldP = vec3(oldP[0], oldP[1], 0); 
			}
			double oldTem = getTemperature(oldP); 
			target.mT (i, j, k) = oldTem; 
			//cout<<"mT"<<target.mT (i, j, k)<<endl; 
		 }
	  } 
   }
   mT = target.mT;    

}

void MACGrid::advectDensity(double dt)
{
   target.mD = mD; 
   for(int k = 0; k < theDim[MACGrid::Z]; k++) 
   {
      for(int j = 0; j < theDim[MACGrid::Y]; j++) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; i++) 
		 {
			vec3 pD = getCenter(i, j, k);    
			vec3 Density = getVelocity(pD); 
			vec3 oldPD = pD - dt*Density;   

			//check boundary for oldPV   
			if(oldPD[0] <= 0)    
			{
				oldPD = vec3(0, oldPD[1], oldPD[2]); 
			}
			if(oldPD[1] <= 0) 
			{
				oldPD = vec3(oldPD[0], 0, oldPD[2]);   
			}
			if(oldPD[2] <= 0) 
			{
				oldPD = vec3(oldPD[0], oldPD[1], 0); 
			}
			double oldDen = getDensity (oldPD);
			target.mD (i, j, k) = oldDen; 
		 }
	  }
   }
   mD = target.mD;     
}

void MACGrid::computeBouyancy(double dt)
{
	// TODO: Calculate bouyancy and store in target.  
	target.mV = mV;
   // Then save the result to our object.    
   double alpha = 5.0; 
   double beta = 15.0; 
   double Tambiant = 20.0; 

   //set boundary for mv           
    for(int k = 0; k < theDim[MACGrid::Z]; k++) 
	{
      for(int j = 1; j < theDim[MACGrid::Y]; j++) 
	  {
         for(int i = 0; i < theDim[MACGrid::X]; i++) 
		 {
		     if((j+1) < theDim[MACGrid::Y])
			 {
				 double Fbouy = -alpha*((mD(i,j,k) + mD(i,j+1,k))/2.0)  + beta*((mT(i,j,k)+ mT(i,j+1, k))/2.0 - Tambiant); 
				 target.mV (i,j,k) = target.mV(i,j,k) + Fbouy*dt; 
			 }
		 }    
	  }  
	}

   mV = target.mV; 
}

void MACGrid::computeVorticityConfinement(double dt)
{
	// TODO: Calculate vorticity confinement forces.
	// Apply the forces to the current velocity and store the result in target.   

    //calculate the memory here 
	int Memlength = theDim[MACGrid::Z]*theDim[MACGrid::Y]*theDim[MACGrid::X];

	vector <vec3> w (Memlength); 
	vector <vec3> wLength (Memlength); 
	vector <vec3> N (Memlength); 
	vector <vec3> Fconf (Memlength); 

	target.mU = mU;
	target.mV = mV;
	target.mW = mW;

	for(int k = 1; k < theDim[MACGrid::Z]; k++) 
	{
      for(int j = 1; j < theDim[MACGrid::Y]; j++) 
	  {
         for(int i = 1; i < theDim[MACGrid::X]; i++) 
		 {
		    int index = i+j*theDim[MACGrid::X]+k*theDim[MACGrid::Y]*theDim[MACGrid::X];
		    w [index] = vec3((target.mW(i, j+1, k)-target.mW(i, j-1, k))/(2*theCellSize)-(target.mV(i, j, k+1)-target.mV(i,j,k-1))/(2*theCellSize), (target.mU(i, j, k+1)-target.mU(i,j,k-1))/(2*theCellSize)-(target.mW(i+1,j,k)-target.mW(i-1,j,k))/(2*theCellSize), (target.mV(i+1, j, k) - target.mV(i-1, j,k))/(2*theCellSize) - (target.mU(i,j+1, k) - target.mU(i,j-1,k))/(2*theCellSize));
		 }    
	  }         
	}  
	
	for(int k = 2; k < theDim[MACGrid::Z]-1; k++)   
	{
      for(int j = 2; j < theDim[MACGrid::Y]-1; j++) 
	  {
         for(int i = 2; i < theDim[MACGrid::X]-1; i++) 
		 {
		    int index = i+j*theDim[MACGrid::X]+k*theDim[MACGrid::Y]*theDim[MACGrid::X];
			int index1 = (i+1) + j*theDim[MACGrid::X] + k*theDim[MACGrid::X]*theDim[MACGrid::Y];   
			int index2 = (i-1) + j*theDim[MACGrid::X] + k*theDim[MACGrid::X]*theDim[MACGrid::Y]; 
			int index3 = i + (j+1)*theDim[MACGrid::X] + k*theDim[MACGrid::X]*theDim[MACGrid::Y]; 
			int index4 = i + (j-1)*theDim[MACGrid::X] + k*theDim[MACGrid::X]*theDim[MACGrid::Y]; 
			int index5 = i + j*theDim[MACGrid::X] + (k+1)*theDim[MACGrid::X]*theDim[MACGrid::Y]; 
			int index6 = i + j*theDim[MACGrid::X] + (k-1)*theDim[MACGrid::X]*theDim[MACGrid::Y]; 

			wLength[index] = vec3((w[index1].Length() - w[index2].Length())/(2*theCellSize), (w[index3].Length() - w[index4].Length())/(2*theCellSize), (w[index5].Length() - w[index6].Length())/(2*theCellSize)); 

		 }    
	  }         
	}    

	
	for(int k = 2; k < theDim[MACGrid::Z]-1; k++) 
	{
      for(int j = 2; j < theDim[MACGrid::Y]-1; j++) 
	  {
         for(int i = 2; i < theDim[MACGrid::X]-1; i++) 
		 {
		    int index = i+j*theDim[MACGrid::X]+k*theDim[MACGrid::Y]*theDim[MACGrid::X];
			if((i>0)&&(j>0)&&(k>0)&&((i+1)<theDim[MACGrid::X])&&((j+1)<theDim[MACGrid::Y])&&((k+1)<theDim[MACGrid::Z]))
			{
				N[index] = wLength[index].Normalize(); 
			}
		 }    
	  }         
	}

	double epislon = 50.00; 
	for(int k = 2; k < theDim[MACGrid::Z]-1; k++) 
	{
      for(int j = 2; j < theDim[MACGrid::Y]-1; j++) 
	  {
         for(int i = 2; i < theDim[MACGrid::X]-1; i++) 
		 {
		    int index = i+j*theDim[MACGrid::X]+k*theDim[MACGrid::Y]*theDim[MACGrid::X];
			if((i>0)&&(j>0)&&(k>0)&&((i+1)<theDim[MACGrid::X])&&((j+1)<theDim[MACGrid::Y])&&((k+1)<theDim[MACGrid::Z]))
			{
				Fconf[index] = epislon*theCellSize*N[index].Cross(w[index]); 
				vec3 temp = Fconf[index]; 

				target.mU (i,j,k) = target.mU(i,j,k) + temp[0]*dt; 
				target.mV (i,j,k) = target.mV(i,j,k) + temp[1]*dt;
				target.mW (i,j,k) = target.mW(i,j,k) + temp[2]*dt;
			} 
		 }    
	  }         
	}  

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
	// Then save the result to our object
    setUpAMatrix(); 

	//initialize d    
	GridData initialNum; 
	initialNum.initialize(0); 
	//GridData& d = initialNum;  
	GridData d; 
	d.initialize(0); 

	for(int k = 0; k < theDim[MACGrid::Z]; k++)
	{
      for(int j = 0; j < theDim[MACGrid::Y]; j++)
	  {
         for(int i = 0; i < theDim[MACGrid::X]; i++) 
		 {
			if((i+1) == theDim[MACGrid::X])
			{
				target.mU ((i+1), j, k) = 0; 
			}
			if((j+1) == theDim[MACGrid::Y])
			{
				target.mV (i, (j+1) , k) = 0; 
			}
			if((k+1) == theDim[MACGrid::Z])
			{
				target.mW (i, j , (k+1)) = 0; 
			}
			if(i == 0)
			{
				target.mU (i, j, k) = 0;
			}
			if(j == 0)
			{
				target.mV (i, j ,k) = 0; 
			}
			if(k == 0)
			{
			    target.mW (i, j , k) = 0; 
			}
		 }
	  }        
   }



    for(int k = 0; k < theDim[MACGrid::Z]; k++)
	{
      for(int j = 0; j < theDim[MACGrid::Y]; j++)
	  {
         for(int i = 0; i < theDim[MACGrid::X]; i++) 
		 {
			 d (i, j, k) = (-1.0)*theCellSize/dt*(target.mU ((i+1), j, k) - target.mU (i, j, k) + target.mV (i, (j+1), k) - target.mV (i, j, k) + target.mW (i, j, (k+1)) - target.mW (i, j, k));    
		 }
	  }        
   }

	//sanity check here 
	bool checkDivergence = false; 

	int i = 240; 
	while(checkDivergence==false)
	{
	
	   checkDivergence = conjugateGradient(AMatrix, target.mP, d, i, 0.01); 
	   i = i+ 20; 
	  // cout<<i<<endl;
	}

   if(checkDivergence)
   {
		for(int k = 0; k < theDim[MACGrid::Z]; k++)       
		{
		  for(int j = 0; j < theDim[MACGrid::Y]; j++)
		  {
			 for(int i = 0; i < theDim[MACGrid::X]; i++)       
			 {
				 if((i+1) < theDim[MACGrid::X])
				 {
					 target.mU (i+1, j, k) = target.mU (i+1, j, k) - dt*(target.mP(i+1, j, k) - target.mP(i, j, k))/(1.0*theCellSize); 
				 }
			 }     
		  }      
		}    
  
		for(int k = 0; k < theDim[MACGrid::Z]; k++)       
		{
		  for(int j = 0; j < theDim[MACGrid::Y]; j++)
		  {
			 for(int i = 0; i < theDim[MACGrid::X]; i++)       
			 {
				 if((j+1) < theDim[MACGrid::Y])
				 {
					 target.mV (i, j+1, k) = target.mV (i, j+1, k) - dt*(target.mP(i, j+1, k) - target.mP(i, j, k))/(1.0*theCellSize); 
				 }
			 }     
		  }      
		} 

		for(int k = 0; k < theDim[MACGrid::Z]; k++)       
		{
		  for(int j = 0; j < theDim[MACGrid::Y]; j++)
		  {
			 for(int i = 0; i < theDim[MACGrid::X]; i++)       
			 {
				 if((k+1) < theDim[MACGrid::Z])
				 {
					 target.mW (i, j, k+1) = target.mV (i, j, k+1) - dt*(target.mP(i, j, k+1) - target.mP(i, j, k))/(1.0*theCellSize); 
				 }
			 }     
		  }     
		} 

		mP = target.mP;
		mU = target.mU;
		mV = target.mV;
		mW = target.mW;
   }else
   {
	    mP = target.mP;
		mU = target.mU;
		mV = target.mV;
		mW = target.mW;
   }
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


void MACGrid::SetParamSigma(double sigma)
{
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
