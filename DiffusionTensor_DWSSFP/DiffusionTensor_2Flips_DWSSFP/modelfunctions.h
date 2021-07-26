#define MACRO template <typename T> __device__ inline

#define gyro 26753.803038

// In DTI model the parameters P are:
// lam1:  0
// lam2:  1
// lam3:  2 
// mu1:   3 
// mu2:   4
// mu3:   5
// theta: 6
// phi:   7
// psi:   8
// S01:   9
// S02:   10

// CFP[0:2] are bvecs 
// CFP[3]   TRs
// CFP[4]   diffGradAmps
// CFP[5]   diffGradDurs
// CFP[6]   flipAngle
// CFP[7]   whichflip  (0=first, 1=second)
// CFP[8]   Noisefloor estimate

// FixP[0] T1
// FixP[1] T2
// FixP[2] B1

MACRO T SSFP_Signal(T adc, T* P, T* CFP, T* FixP, T E1, T E2, T sa, T ca)
{
  T S0;
  if(CFP[7]==0){
      S0 = P[9];
  }
  else{
      S0 = P[10];
  }

  T qval = (T)gyro*CFP[4]*CFP[5]/(T)10.0;
  // Iso
  T A1    = exp_gpu( -qval*qval*CFP[3]*adc );
  T A2    = exp_gpu( -qval*qval*CFP[5]*adc );
  T A2_03 = exp_gpu( -qval*qval*CFP[5]*adc/(T)3.0 );

  T s     = E2*A1/A2_03/A2_03/A2_03/A2_03*((T)1.0-E1*ca)+E2/A2_03*(ca-E1);
  T r     = (T)1.0 - E1*ca+E2*E2*A1*A2_03*(ca-E1);
  T K     = ((T)1.0-E1*A1*ca-E2*E2*A1*A1/A2_03/A2_03*(E1*A1-ca))/(E2*A1/A2_03/A2_03/A2_03/A2_03*((T)1.0+ca)*((T)1.0-E1*A1));

  T F1    = K - sqrt_gpu(K*K-A2*A2);
  T Mminus_top = -((T)1.0-E1)*E2/A2_03/A2_03*(F1-E2*A1*A2_03*A2_03)*sa;
  T Mminus_bottom=r-F1*s;
  T signal = sqrt_gpu(fabs_gpu(S0*Mminus_top/Mminus_bottom)*fabs_gpu(S0*Mminus_top/Mminus_bottom)+CFP[8]*CFP[8]);

  return signal;
}

MACRO T Compute_ADC(T* P, T* CFP)
{  


  T lam1  = P[0];
  T lam2  = P[1];
  T lam3  = P[2];
  T mu1   = P[3];
  T mu2   = P[4];
  T mu3   = P[5];
  T theta = P[6];
  T phi   = P[7];
  T psi   = P[8];

  T v1x,v1y,v1z,v2x,v2y,v2z,v3x,v3y,v3z;
  T s = sin_gpu(theta);

  v1x = s*cos_gpu(phi);
  v1y = s*sin_gpu(phi);
  v1z = cos_gpu(theta);

  if( (fabs_gpu(s)>0) ){
    v2x = (-cos_gpu(psi)*v1y - sin_gpu(psi)*v1x*v1z)/sqrt_gpu(s*s);
    v2y = (cos_gpu(psi)*v1x - sin_gpu(psi)*v1y*v1z)/sqrt_gpu(s*s);
    v2z = sin_gpu(psi)*sqrt(s*s);
  }
  else{
    v2x = -cos_gpu(psi);
    v2y = 0.0;
    v2z = sin_gpu(psi);
  }
  v3x = v1y*v2z - v2y*v1z;
  v3y = v2x*v1z - v1x*v2z;
  v3z = v1x*v2y - v2x*v1y;
  
  T Dxx,Dyy,Dzz,Dxy,Dxz,Dyz;
  if(CFP[7]==0){
    Dxx = lam1*v1x*v1x + lam2*v2x*v2x + lam3*v3x*v3x;
    Dxy = lam1*v1x*v1y + lam2*v2x*v2y + lam3*v3x*v3y;
    Dxz = lam1*v1x*v1z + lam2*v2x*v2z + lam3*v3x*v3z;
    Dyy = lam1*v1y*v1y + lam2*v2y*v2y + lam3*v3y*v3y;
    Dyz = lam1*v1y*v1z + lam2*v2y*v2z + lam3*v3y*v3z;
    Dzz = lam1*v1z*v1z + lam2*v2z*v2z + lam3*v3z*v3z;
  }
  else{
    Dxx = mu1*v1x*v1x + mu2*v2x*v2x + mu3*v3x*v3x;
    Dxy = mu1*v1x*v1y + mu2*v2x*v2y + mu3*v3x*v3y;
    Dxz = mu1*v1x*v1z + mu2*v2x*v2z + mu3*v3x*v3z;
    Dyy = mu1*v1y*v1y + mu2*v2y*v2y + mu3*v3y*v3y;
    Dyz = mu1*v1y*v1z + mu2*v2y*v2z + mu3*v3y*v3z;
    Dzz = mu1*v1z*v1z + mu2*v2z*v2z + mu3*v3z*v3z;
  }

  T adc =      CFP[0]*CFP[0]*Dxx
	   +   CFP[1]*CFP[1]*Dyy 
           +   CFP[2]*CFP[2]*Dzz
           + 2*CFP[0]*CFP[1]*Dxy
           + 2*CFP[0]*CFP[2]*Dxz
           + 2*CFP[1]*CFP[2]*Dyz ;

  return adc;
}

MACRO T Predicted_Signal(
			 int npar, // Number of Parameters to estimate
			 T* P, // Estimated parameters
			 T* CFP, // Fixed Parameters common to all the voxels
			 T* FixP) // Fixed Parameters for each voxel
{

  // Relaxation terms
  T E1 = exp_gpu(-CFP[3]/ ((T)1e-3*FixP[0]));
  T E2 = exp_gpu(-CFP[3]/ ((T)1e-3*FixP[1]));

  // Flipping
  T sa = sin_gpu(CFP[6]*FixP[2]*M_PI/(T)180.0);
  T ca = cos_gpu(CFP[6]*FixP[2]*M_PI/(T)180.0);

  // Diffusion terms  
  T adc = Compute_ADC(P,CFP);
  T pred_signal = SSFP_Signal(adc,P,CFP,FixP,E1,E2,sa,ca);

  return pred_signal;
}

// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
			   int npar, // Number of Parameters to estimate
			   T* P) // Estimated parameters
{
  
  // preserve eigenvalue ordering to ensure eigenvectors are shared across flip angles
  if( P[0]<P[1] || P[1]<P[2] || P[0]<P[2] || P[3]<P[4] || P[4]<P[5] || P[3]<P[5] ){ 
    return false;
  }
  return true;
}

// Partial derivatives respect each model parameter
MACRO void Partial_Derivatives(
			       int npar, // Number of Parameters to estimate
			       T* P, // Estimated parameters, use P*
			       T* CFP, // Fixed Parameters common to all the voxels
			       T* FixP, // Fixed Parameters for each voxel
			       T* derivatives) // Derivative respect each model estimated parameter
{
  
}

// Constraints run after LevenbergMarquardt (if LevenbergMarquardt is used)
MACRO void FixConstraintsLM(	
			    int npar, // Number of Parameters to estimate
			    T* P) // Estimated parameters
{


}

MACRO T custom_priors(
			int id_p,   // the number of parameter in the model (starts at 0)
			T* P, 		// Estimated parameters
	int nmeas, // Number of measurements per voxel
			 T* CFP, 	// Fixed Parameters common to all the voxels for all measurements !!
			 T* FixP) 	// Fixed Parameters for each voxel
{
	return -log_gpu(P[id_p]);
}

