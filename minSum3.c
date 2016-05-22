#include "mex.h"
#include <math.h>

#define INF 1e20
#define printf mexPrintf

void minSum3(double *cHat, double *gamma,double *H, int m, int n, int maxIterations);
void minSumOperation(double *y, double *z, int n, int *d, int nd);

/* mex Function */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] ) {

/* input : vReceived = prhs[0], received vector
 *         H = prhs[1], H transpose matrix
 *         maxIterations = prhs[2], number of maxIterations 
 *
 * output: vHat = plhs[0], hard decision output
 */         
   
  	//if(nrhs != 3) mexErrMsgTxt("Required 3 arguments as follow: gamma, H, maxIterations");
   
    /* Create pointer to the inputs */
    double *gamma 			= mxGetPr(prhs[0]);      
    double *H 				= mxGetPr(prhs[1]);
    double maxIterations 	= mxGetScalar(prhs[2]);
   
    //get number of rows
    int m = mxGetM(prhs[1]);

    //get number of columns 
    int n = mxGetN(prhs[1]);

    //set the output pointer to the output vector 
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);   
   
    //create a pointer to a copy ot the output vector cHat 
    double *cHat = mxGetPr(plhs[0]);
   
    //calling function
    minSum3(cHat, gamma, H, m, n, (int) maxIterations);
}

void minSum3(double *cHat, double *gamma, double *H, int m, int n, int maxIterations) {

	//p pointer to the sparse matrix
	//init with m*n as number of ones (Worst case)
	int *p = (int *) calloc(m*n, sizeof(int)); 
	
	//d collects number of checks per row of H
	int d[m];
	
	//construction of P. This may take a while
	int nnz = 0; //when cycle finish in nnz i have number of nonzeros
	for (int i = 0; i < m;i++) {
		d[i] = 0; //init d[i] value
		for (int j= 0; j < n;j++) {
			if (*(H+j*m+i) != 0) {
				*(p+nnz) = j;
				d[i]++;
				nnz++;
			}
		}
	}
	
	//reallocating memory since we know real nnz
	p = (int *) realloc(p, nnz*sizeof(int));
	
	//init variable node value from LogLikelihoodRatio gamma
	double z[nnz];
	
	for (int k = 0; k < nnz;k++) {
		int curr = *(p+k);
		z[k] = *(gamma+curr);
	}
	
	//init variables
	double y[nnz];
	double u[n];
	
	for (int iter = 0; iter<maxIterations;iter++) {
	
		//updating check nodes
		minSumOperation(&y[0], &z[0], nnz, &d[0], m);
		
		//init decision vector with variable nodes
		double u[n];
		for (int k = 0; k < n; k++) {
			u[k] = *(gamma+k);
		}
		
		//updating decision vector with check nodes
		for (int k = 0; k < nnz; k++) {
			int curr = *(p+k);
			u[curr] += y[k];
		}
		
		//decision
		for (int k= 0; k< n;k++) {
			if (u[k] < 0) {
				*(cHat+k) = 1;
			} else {
				*(cHat+k) = 0;
			}
		}
	
		//update variable nodes
		for (int k = 0; k < nnz;k++) {
			int curr = *(p+k);
			z[k] = u[curr]-y[k];
		}
	}
	
	//free memory
	free(p);
	
}

void minSumOperation(double *y, double *z, int n, int *d, int nd) {
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