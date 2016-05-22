#include "mex.h"
#include <math.h>

#define INF 1e20
#define printf mexPrintf

void minSumOperation(double *y, double *z, int n, double *d, int nd);

/* mex Function */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] ) {

/* input : vReceived = prhs[0], received vector
 *         H = prhs[1], H transpose matrix
 *         maxIterations = prhs[2], number of maxIterations 
 *
 * output: vHat = plhs[0], hard decision output
 */         

    double *z, *y;
    int n, nd;
   
  	if(nrhs != 2) mexErrMsgTxt("Two argument required");
   
    /* Create pointer to the inputs */
    z = mxGetPr(prhs[0]);
    
    double *d = mxGetPr(prhs[1]);
    
      /* Get the # of column of the H transpose */
    n = mxGetM(prhs[0]);
    
    nd = mxGetM(prhs[1]);
    
    double count = 0;
  	for (int i=0; i < nd; i++) {
  		count += *(d+i);
  	}
  	
  	if (count != n) mexErrMsgTxt("Incompatible vectors");

    /* Set the output pointer to the output vector */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);   
   
    /* Create a pointer to a copy ot the output vector LOut */
    y = mxGetPr(plhs[0]);
   
    /* Call the function */
    minSumOperation(y, z, n, d, nd);
}

/* minSum function */
void minSumOperation(double *y, double *z, int n, double *d, int nd) {
	int sign;
	int k = 0;
	
	for (int i = 0; i < nd; i++) {
		int length = (int) *(d+i);
		
		for (int j = 0; j < length; j++) {
			double mintmp = INF;
			sign = 1;
		
			for (int yi = 0; yi < length; yi++) {
				if (yi!=j) {
					int signtmp = (*(z+yi+k))<0?-1:1;
					sign *= signtmp;
					
					if (fabs(*(z+yi+k))<mintmp) {
						mintmp = fabs(*(z+yi+k));
					}
				}
			}
			*(y+k+j) = mintmp*sign;
		}
		
		 k+= length;
	}
}