/****************************************************************************\
*									     *
*			       Vector addition				     *
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
void fvadd(double* o, const double* s, const double* t, int n);
void fvadds(double* o, const double* s, const double* t, int n);
void fvmov(double* o, const double* s, int n);
void fvmovs(double* o, const double* s, int n);


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{	int m, n;
	if (nlhs > 1) mexErrMsgTxt("Too many output arguments.");
	switch(nrhs)
	{default:
		mexErrMsgTxt("Incorrect number of arguments.");
	 case 0:
		mexPrintf("\nVector addition.\n\n\tMarcel Leutenegger © 4.12.2005\n\n");
		break;
	 case 2:
		m=mxGetN(S);
		n=mxGetN(T);
		O=mxCreateDoubleMatrix(0,0,mxREAL);
		if (m > 0 && n > 0)
		{	double* o;
			const int* d=mxGetDimensions(S);
			const int* e=mxGetDimensions(T);
			if (*d == 3 && *e == 3)
			{	mxComplexity cs=mxIsComplex(S);
				mxComplexity ct=mxIsComplex(T);
				if (m == 1)	// vadd(vector,matrix)
				{	mxSetDimensions(O,e,mxGetNumberOfDimensions(T));
					mxSetPr(O,o=mxMalloc(m=n*3*sizeof(double)));
					fvadds(o,mxGetPr(S),mxGetPr(T),n);
					if (cs || ct)
					{	mxSetPi(O,o=mxMalloc(m));
						if (cs)
						{	if (ct)
								fvadds(o,mxGetPi(S),mxGetPi(T),n);
							else	fvmovs(o,mxGetPi(S),n);
						}
						else	fvmov(o,mxGetPi(T),n);
					}
					break;
				}
				if (n == 1)	// vadd(matrix,vector)
				{	mxSetDimensions(O,d,mxGetNumberOfDimensions(S));
					mxSetPr(O,o=mxMalloc(n=m*3*sizeof(double)));
					fvadds(o,mxGetPr(T),mxGetPr(S),m);
					if (cs || ct)
					{	mxSetPi(O,o=mxMalloc(n));
						if (cs)
						{	if (ct)
								fvadds(o,mxGetPi(T),mxGetPi(S),m);
							else	fvmov(o,mxGetPi(S),m);
						}
						else	fvmovs(o,mxGetPi(T),m);
					}
					break;
				}
				if (m == n)	// vadd(matrix,matrix)
				{	m=mxGetNumberOfDimensions(S);
					if (fdim(m,mxGetNumberOfDimensions(T),d,e))
					{	mxSetDimensions(O,d,m);
						mxSetPr(O,o=mxMalloc(m=n*3*sizeof(double)));
						fvadd(o,mxGetPr(S),mxGetPr(T),n);
						if (cs || ct)
						{	mxSetPi(O,o=mxMalloc(m));
							if (cs)
							{	if (ct)
									fvadd(o,mxGetPi(S),mxGetPi(T),n);
								else	fvmov(o,mxGetPi(S),n);
							}
							else	fvmov(o,mxGetPi(T),n);
						}
						break;
					}
				}
			}
			mexErrMsgTxt("Incompatible dimensions.");
		}
	}
}
