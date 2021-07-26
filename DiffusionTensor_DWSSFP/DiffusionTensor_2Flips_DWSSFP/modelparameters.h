#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

typedef double MyType;

///// Edit this /////
#define NPARAMS 11
// bvecs,TRs,diffGradAmps,diffGradDurs,flipangle,flipmask,noisefloor
#define NCFP 7
// T1,T2,B1
#define NFIXP 3
/////////////////////

///// Do not edit this /////
struct MODEL
{
	static int CFP_size[NCFP];
	static int FixP_size[NFIXP];
};
////////////////////////////

#endif

