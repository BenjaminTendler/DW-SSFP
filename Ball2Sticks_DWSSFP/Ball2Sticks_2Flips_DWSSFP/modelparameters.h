#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

// Do not forget quotation marks for specifing the model files.
// Use absolute paths.

typedef double MyType;

// d1,d2,f1,th1,ph1,f2,th2,ph2,S01,S02 
#define NPARAMS 10
// bvecs,bvals,TRs,diffGradAmps,diffGradDurs,flipangle,flipMask,noisefloor
#define NCFP 8
// T1,T2,B1
#define NFIXP 3

struct MODEL 
{
	static int CFP_size[NCFP];
	static int FixP_size[NFIXP];
};

#endif

