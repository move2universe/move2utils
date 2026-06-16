#include <R.h>
#include <math.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
SEXP llBGBvar(SEXP sigm, SEXP paraOrth)
{
	SEXP rDim, ans;
	int I;
//		warning("%f",sigma);
	rDim=getAttrib(sigm, R_DimSymbol);
	I=INTEGER(rDim)[0];
	PROTECT(ans=allocVector(REALSXP,1));
	REAL(ans)[0]=0;
//		warning("%f",REAL(ans)[0]);
//		printf("asdf %f", REAL(ans)[0]);
	for(int i =0; i<I;i++)
	{

		REAL(ans)[0]+=log(1/(2*3.141593*sqrt(REAL(sigm)[i])*sqrt(REAL(sigm)[i+I])))-0.5 * ((pow(REAL(paraOrth)[i],2)/REAL(sigm)[i]) + (pow(REAL(paraOrth)[i+I],2)/REAL(sigm)[i+I]));
	}
	UNPROTECT(1);
	return(ans);
}
