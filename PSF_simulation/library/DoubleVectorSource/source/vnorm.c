/****************************************************************************\
*									     *
*				  Vector norm				     *
*									     *
*									     *
*			Marcel Leutenegger © 4.12.2005			     *
*									     *
\****************************************************************************/

#include "mex.h"
#define	O	plhs[0]
#define	S	prhs[0]

void fvnorm(double* or, const double* sr, const double* si, int n);


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{	int n;
	double* or;
	const int* d;
	if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
	switch(nrhs)
	{default:
		mexErrMsgTxt("Incorrect number of arguments.");
	 case 0:
		mexPrintf("\nVector norm.\n\n\tMarcel Leutenegger © 4.12.2005\n\n");
		break;
	 case 1:
		d=mxGetDimensions(S);
		if (*d != 3) mexErrMsgTxt("Incompatible dimensions.");
		O=mxCreateDoubleMatrix(0,0,mxREAL);
		n=mxGetN(S);
		if (n > 0)
		{	mxSetDimensions(O,d,mxGetNumberOfDimensions(S));
			mxSetPr(O,or=mxMalloc(n*sizeof(double)));
			mxSetM(O,1);
			fvnorm(or,mxGetPr(S),mxGetPi(S),n);
		}
	}
}
