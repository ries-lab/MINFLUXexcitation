/****************************************************************************\
*									     *
*			      Vector dot product			     *
*									     *
*									     *
*			Marcel Leutenegger © 4.12.2005			     *
*									     *
\****************************************************************************/

#include "mex.h"
#define	O	plhs[0]
#define	S	prhs[0]
#define	T	prhs[1]

bool fdim(int m, int n, const int* d, const int* e);
void fvdot(double* or, double* oi, const double* sr, const double* si, const double* tr, const double* ti, int n);
void fvdots(double* or, double* oi, const double* sr, const double* si, const double* tr, const double* ti, int n);
void fvdott(double* or, double* oi, const double* sr, const double* si, const double* tr, const double* ti, int n);


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{	int m, n;
	if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
	switch(nrhs)
	{default:
		mexErrMsgTxt("Incorrect number of arguments.");
	 case 0:
		mexPrintf("\nVector dot product.\n\n\tMarcel Leutenegger © 4.12.2005\n\n");
		break;
	 case 2:
		m=mxGetN(S);
		n=mxGetN(T);
		O=mxCreateDoubleMatrix(0,0,mxREAL);
		if (m > 0 && n > 0)
		{	double* or;
			double* oi=NULL;
			const int* d=mxGetDimensions(S);
			const int* e=mxGetDimensions(T);
			if (*d == 3 && *e == 3)
			{	mxComplexity cplx=mxIsComplex(S) || mxIsComplex(T);
				if (m == 1)	// vdot(vector,matrix)
				{	mxSetDimensions(O,e,mxGetNumberOfDimensions(T));
					mxSetPr(O,or=mxMalloc(m=n*sizeof(double)));
					if (cplx) mxSetPi(O,oi=mxMalloc(m));
					mxSetM(O,1);
					fvdots(or,oi,mxGetPr(S),mxGetPi(S),mxGetPr(T),mxGetPi(T),n);
					break;
				}
				if (n == 1)	// vdot(matrix,vector)
				{	mxSetDimensions(O,d,mxGetNumberOfDimensions(S));
					mxSetPr(O,or=mxMalloc(n=m*sizeof(double)));
					if (cplx) mxSetPi(O,oi=mxMalloc(n));
					mxSetM(O,1);
					fvdott(or,oi,mxGetPr(S),mxGetPi(S),mxGetPr(T),mxGetPi(T),m);
					break;
				}
				if (m == n)	// vdot(matrix,matrix)
				{	m=mxGetNumberOfDimensions(S);
					if (fdim(m,mxGetNumberOfDimensions(T),d,e))
					{	mxSetDimensions(O,d,m);
						mxSetPr(O,or=mxMalloc(m=n*sizeof(double)));
						if (cplx) mxSetPi(O,oi=mxMalloc(m));
						mxSetM(O,1);
						fvdot(or,oi,mxGetPr(S),mxGetPi(S),mxGetPr(T),mxGetPi(T),n);
						break;
					}
				}
			}
			mexErrMsgTxt("Incompatible dimensions.");
		}
	}
}
