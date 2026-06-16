#include <R.h>
#include <math.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

SEXP dbbmm2(SEXP x, SEXP y, SEXP s, SEXP t, SEXP locEr,SEXP xGrid, SEXP yGrid,SEXP tStep, SEXP ext2, SEXP interest)
{
	double *rans, xRes, yRes, x0, y0,ZTZ,  ext, *xx, *xy,*xs, *xxGrid, *xyGrid, alpha, dT,sigma, ti, *xt, mux, muy, *xlocEr;  
	int xEnd, xStart, yEnd, yStart;
	R_len_t i,j,k,nLoc=length(x), nXGrid=length(xGrid), nYGrid=length(yGrid);
	PROTECT(ext2 = coerceVector(ext2, REALSXP));
	ext=REAL(ext2)[0];
	PROTECT(interest= coerceVector(interest, LGLSXP));
	PROTECT(xGrid = coerceVector(xGrid, REALSXP));
	PROTECT(yGrid = coerceVector(yGrid, REALSXP));
	PROTECT(x = coerceVector(x, REALSXP));
	PROTECT(y = coerceVector(y, REALSXP));
	PROTECT(s = coerceVector(s, REALSXP));
	PROTECT(t = coerceVector(t, REALSXP));
	PROTECT(tStep = coerceVector(tStep, REALSXP));
	PROTECT(locEr = coerceVector(locEr, REALSXP));
	dT=REAL(tStep)[0];
	
	x0 = REAL(xGrid)[0];
	y0 = REAL(yGrid)[0];
	xt = REAL(t);
	xlocEr=REAL(locEr);
	xx=REAL(x);
	xy=REAL(y);
	xyGrid=REAL(yGrid);
	xxGrid=REAL(xGrid);
	xRes = xxGrid[1]-xxGrid[0];
	yRes = xyGrid[1]-xyGrid[0];
	xs=REAL(s);
	
//	warning("%f", ti);
SEXP ans;
	PROTECT(ans=allocMatrix(REALSXP, nYGrid, nXGrid));
	rans = REAL(ans);
	for(i = 0; i < nXGrid; i++) {
		for(j = 0; j < nYGrid; j++)
			rans[i* nYGrid+ j] = 0.0;
	}
//	for(k=0; k<nLoc;k++)
	k=0; //k counts the segment were in
	for(ti=xt[0]+(fmod((xt[nLoc-1]-xt[0]), dT)/2.0); ti<=xt[nLoc-1]; ti=ti+dT)//check mod seems not to work for exact things
	{
		while(xt[k+1]<ti)
		{
			k++;
		}
		if(LOGICAL(interest)[k]){// only calculate for segments of interest
  		alpha=(ti-xt[k])/(xt[k+1]-xt[k]);
  		mux=xx[k]+(xx[k+1]-xx[k])*alpha;
	  	muy=xy[k]+(xy[k+1]-xy[k])*alpha;
  		sigma=(xt[k+1]-xt[k])*alpha*(1-alpha)*xs[k]+pow(1-alpha,2)*pow(xlocEr[k],2)+pow(alpha,2)*pow(xlocEr[k+1],2);

  		xStart=(int)floor(((mux-x0)/xRes)-((sqrt(sigma)*ext)/xRes));
	  	xEnd=(int)ceil(((mux-x0)/xRes)+((sqrt(sigma)*ext)/xRes));
		  yEnd=(int)nYGrid-floor(((muy-y0)/yRes)-((sqrt(sigma)*ext)/yRes));
  		yStart=(int)nYGrid-ceil(((muy-y0)/yRes)+((sqrt(sigma)*ext)/yRes));
	  	if(xStart<0){
//			error("Lower x grid not large enough %i %.16f s=%f dt=%f a=%f k=%i l %f %f t %f %f, mux=%f, X0=%f, xRes=%f, ext=%f, var=%f", xStart, ti, sigma,  xt[k+1]-xt[k], alpha, k, xlocEr[k], xlocEr[k+1], xt[k], xt[k+1], mux, x0, xRes, ext, xs[k]);
		  	error("Lower x grid not large enough, consider extending the raster in that direction or enlarging the ext argument");
	  	}
	  	if(xEnd>nXGrid){
	  	  error("Higher x grid not large enough, consider extending the raster in that direction or enlarging the ext argument");
		  }
  		if(yEnd>nYGrid){
	  		error("Lower y grid not large enough, consider extending the raster in that direction or enlarging the ext argument");
		  }
	  	if(yStart<0){
	  		error("Higher y grid not large enough, consider extending the raster in that direction or enlarging the ext argument");
		  }
	  	R_CheckUserInterrupt();
//``warning("my %f %f %f %f", muy, xyGrid[yStart], xyGrid[yEnd], sigma);
	//	tmp=0;
	//	for(i = xStart ; i<xEnd; i++){
	//		for(j = yStart ; j<yEnd; j++){
	//			ZTZ=pow(xxGrid[i]-mux,2)+pow(xyGrid[nYGrid-j-1]-muy,2);//check nYGrid for 0/1 starting errors
	//			tmp=tmp+((1/(2*3.141593*sigma))*exp(-ZTZ/(2*sigma)));
	//		}}
		  for(i = xStart ; i<=xEnd; i++){
		  	for(j = yStart ; j<=yEnd; j++){
			  	ZTZ=pow(xxGrid[i]-mux,2)+pow(xyGrid[nYGrid-j-1]-muy,2);//check nYGrid for 0/1 starting errors
				  rans[i*nYGrid +j ]=rans[i*nYGrid +j ]+((1/(2*3.141593*sigma))*exp(-ZTZ/(2*sigma)))*xRes*yRes;
		  	}
		  }
	  }
	}
	UNPROTECT(11);
	return(ans);
}
SEXP bgb(SEXP x, SEXP y, SEXP sPara,SEXP sOrth, SEXP t, SEXP locEr,SEXP xGrid, SEXP yGrid, SEXP dTT, SEXP ext2)
{
	double *rans, xRes, yRes, x0, y0,  ext, *xx, *xy,*xsPara, *xsOrth,tmp, *xxGrid, *xyGrid, alpha, sigmaOrth, sigmaPara, ti, *xt, mux, muy, *xlocEr, deltaPara,deltaOrth,A,B,C, dT;  
	int xEnd, xStart, yEnd, yStart;
	PROTECT(ext2 = coerceVector(ext2, REALSXP));
	ext=REAL(ext2)[0];
	PROTECT(dTT=coerceVector(dTT, REALSXP));
	dT=REAL(dTT)[0];
//	PROTECT(maxInt=coerceVector(maxInt, REALSXP));
//	maxT=REAL(maxInt)[0];
	R_len_t i,j,k,nLoc=length(x), nXGrid=length(xGrid), nYGrid=length(yGrid);
	PROTECT(xGrid = coerceVector(xGrid, REALSXP));
	PROTECT(yGrid = coerceVector(yGrid, REALSXP));
	PROTECT(x = coerceVector(x, REALSXP));
	PROTECT(y = coerceVector(y, REALSXP));
	PROTECT(sPara = coerceVector(sPara, REALSXP));
	PROTECT(sOrth = coerceVector(sOrth, REALSXP));
	PROTECT(t = coerceVector(t, REALSXP));
	PROTECT(locEr = coerceVector(locEr, REALSXP));

	x0 = REAL(xGrid)[0];
	y0 = REAL(yGrid)[0];
	xt = REAL(t);
	xlocEr=REAL(locEr);
	xx=REAL(x);
	xy=REAL(y);
	xyGrid=REAL(yGrid);
	xxGrid=REAL(xGrid);
	xRes = xxGrid[1]-xxGrid[0];
	yRes = xyGrid[1]-xyGrid[0];
	xsPara=REAL(sPara);
	xsOrth=REAL(sOrth);
	SEXP ans;
	PROTECT(ans=allocMatrix(REALSXP, nYGrid, nXGrid));
	rans = REAL(ans);
	for(i = 0; i < nXGrid; i++) {
		for(j = 0; j < nYGrid; j++)
			rans[i* nYGrid+ j] = 0.0;
	}
//	for(k=0; k<nLoc;k++)
	k=0; //k counts the segment were in
//	warning("%f",fmod((xt[nLoc-1]-xt[0]), dT));
//	warning("%f",xt[0]+fmod((xt[nLoc-1]-xt[0]), dT));
//	warning("%f",xt[0]);
//	warning("%f",xt[1]);
	for(ti=xt[0]+(fmod((xt[nLoc-1]-xt[0]), dT)/2.0); ti<=xt[nLoc-1]; ti=ti+dT)//check mod seems not to work for exact things
	{
		while(xt[k+1]<ti)
		{
			k++;
		}
//		if((xt[k+1]-xt[k])<maxT){
		R_CheckUserInterrupt();
		alpha=(ti-xt[k])/(xt[k+1]-xt[k]);
		mux=xx[k]+(xx[k+1]-xx[k])*alpha;
		muy=xy[k]+(xy[k+1]-xy[k])*alpha;
		sigmaPara=sqrt((xt[k+1]-xt[k])*alpha*(1-alpha)*pow(xsPara[k],2)+pow(1-alpha,2)*pow(xlocEr[k],2)+pow(alpha,2)*pow(xlocEr[k+1],2));
		sigmaOrth=sqrt((xt[k+1]-xt[k])*alpha*(1-alpha)*pow(xsOrth[k],2)+pow(1-alpha,2)*pow(xlocEr[k],2)+pow(alpha,2)*pow(xlocEr[k+1],2));

		/* calculate over which part of the grid to calculate */
		xStart=(int)floor(((mux-x0)/xRes)-(((fmax2(sigmaPara, sigmaOrth))*ext)/xRes));
		xEnd=(int)ceil(((mux-x0)/xRes)+(((fmax2(sigmaPara, sigmaOrth))*ext)/xRes));
		yEnd=(int)nYGrid-floor(((muy-y0)/yRes)-(((fmax2(sigmaPara, sigmaOrth))*ext)/yRes));
		yStart=(int)nYGrid-ceil(((muy-y0)/yRes)+(((fmax2(sigmaPara, sigmaOrth))*ext)/yRes));
		/* check if the grid is large enough */
		if(xStart<0){
//			error("Lower x grid not large enough %i %.16f spara=%f sorth=%f dt=%f a=%f k=%i l %f %f t %f %f, mux=%f, X0=%f, xRes=%f, ext=%f, var=%f %f", xStart, ti, sigmaPara, sigmaOrth,  xt[k+1]-xt[k], alpha, k, xlocEr[k], xlocEr[k+1], xt[k], xt[k+1], mux, x0, xRes, ext, xsPara[k], xsOrth[k]);
			error("The raster does not extent far enough in the X dimension towards the left side");
		}
		if(xEnd>nXGrid){
			error("The raster does not extent far enough in the X dimension towards the right side");
		}
		if(yEnd>nYGrid){
			error("The raster does not extent far enough in the Y dimension towards the lower side");
		}
		if(yStart<0){
			error("The raster does not extent far enough in the Y dimension towards the upper side");
		}
//		/* loop over grid */
		for(i = xStart ; i<=xEnd; i++){
			for(j = yStart ; j<=yEnd; j++){
				C=sqrt(pow(mux-xx[k+1],2)+pow(muy-xy[k+1],2));
				B=sqrt(pow(xxGrid[i]-xx[k+1],2)+pow(xyGrid[nYGrid-j-1]-xy[k+1],2));
				A=sqrt(pow(mux-xxGrid[i],2)+pow(muy-xyGrid[nYGrid-j-1],2));
				tmp=(B*B-A*A-C*C)/(2*A*C);
				if(tmp>1)
				{ tmp=1;
				}
				if(tmp< -1)
				{ tmp=-1;
				}
				if(C==0){
					tmp=sqrt(.5);
				}
				deltaOrth=A*sin(acos(tmp));
				deltaPara=A*cos(acos(tmp));
				if(A==0){
					deltaPara=0;
					deltaOrth=0;
				}
//				warning("a %.20f b %f c %f dp %f do %f", A, B,C, deltaPara,deltaOrth);
//				warning("a %.20f b %f c %f dp %f do %f", A, B,C, deltaPara,deltaOrth);
//				warning("%f",cos(acos((B*B-A*A-C*C)/(2*A*C))));
//				warning("%f",(acos((B*B-A*A-C*C)/(2*A*C))));
//				warning("%f",(((B*B-A*A-C*C)/(2*A*C))));
				//if(ISNAN(deltaPara)){
				//warning("%f",tmp);
				//}
				rans[i*nYGrid +j ]=rans[i*nYGrid +j ]+((1/(2*3.141593*(sigmaPara)*(sigmaOrth)))*exp(-0.5*((pow(deltaOrth,2)/pow(sigmaOrth,2))+(pow(deltaPara,2)/pow(sigmaPara,2))))*xRes*yRes);
//				warning("%f", rans[i*nYGrid+j]);
//				warning("a %f",((1/(2*3.14*(sigmaPara)*(sigmaOrth)))*exp(-0.5*((pow(deltaOrth,2)/pow(sigmaOrth,2))+(pow(deltaPara,2)/pow(sigmaPara,2))))*xRes*yRes));
//				rans[i*nYGrid+j]
	}}}
	UNPROTECT(11);
	return(ans);
}
