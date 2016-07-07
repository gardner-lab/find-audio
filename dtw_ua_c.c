/**
 * Copyright (C) 2013 Quan Wang <wangq10@rpi.edu>,
 * Signal Analysis and Machine Perception Laboratory,
 * Department of Electrical, Computer, and Systems Engineering,
 * Rensselaer Polytechnic Institute, Troy, NY 12180, USA
 *
 * Modified by L. Nathan Perkins to support unannchored `s`.
 */

/** 
 * This is the C/MEX code of dynamic time warping of two signals
 *
 * compile: 
 *     mex dtw_ua_c.c
 *
 * usage:
 *     d=dtw_ua_c(s,t)
 *     where s is unanchored signal 1, t is signal 2
 */

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
//#if __APPLE__
//#include <Accelerate/Accelerate.h>
//#endif

double vectorDistance(double *s, double *t, int ns, int nt, int k, int i, int j) {
    double result=0;
//#if __APPLE__
//    // surprisingly, this appears to be slower...
//    vDSP_distancesqD(s + k * i, 1, t + k * j, 1, &result, k);
//#else
    double ss, tt;
    int x;
    for (x = 0; x < k; x++) {
        ss = s[k * i + x];
        tt = t[k * j + x];
        result += ((ss - tt) * (ss - tt));
    }
//#endif
    result = sqrt(result);
    return result;
}

void dtw_ua_c(double *s, double *t, int ns, int nt, int k, double *dp) {
    double *D[2];
    int i, j;
    int j1, j2;
    double cost, temp;
    
    // create D
    D[0] = (double *)mxMalloc((nt + 1) * sizeof(double));
    D[1] = (double *)mxMalloc((nt + 1) * sizeof(double));
    
    // initialization
    for (j = 0; j <= nt; j++) {
        D[0][j] = 0;
    }
    
    // dynamic programming
    for (i = 1; i <= ns; i++) {
        D[i % 2][0] = DBL_MAX;
        
        for (j = 1; j <= nt; j++) {
            // calculate norm
            cost = vectorDistance(s, t, ns, nt, k, i - 1, j - 1);
            
            temp = D[(i - 1) % 2][j];
            if (D[i % 2][j - 1] < temp) temp = D[i % 2][j - 1];
            if(D[(i - 1) % 2][j - 1] < temp) temp = D[(i - 1) % 2][j - 1];
            
            D[i % 2][j] = cost + temp;
        }
    }
    
    /* copy to output vector */
    memcpy(dp, D[ns % 2] + 1, nt * sizeof(double));
    
    // free D
    mxFree(D[0]);
    mxFree(D[1]);
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *s,*t;
    int w;
    int ns,nt,k;
    double *dp;
    double *result;
    
    /*  check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:dtw_ua_c:invalidNumInputs", "Two inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt( "MATLAB:dtw_ua_c:invalidNumOutputs", "dtw_ua_c: One output required.");
    }
    
    /*  create a pointer to the input matrix s */
    s = mxGetPr(prhs[0]);
    
    /*  create a pointer to the input matrix t */
    t = mxGetPr(prhs[1]);
    
    /*  get the dimensions of the matrix input s */
    ns = mxGetN(prhs[0]);
    k = mxGetM(prhs[0]);
    
    /*  get the dimensions of the matrix input t */
    nt = mxGetN(prhs[1]);
    if(mxGetM(prhs[1]) != k) {
        mexErrMsgIdAndTxt("MATLAB:dtw_ua_c:dimNotMatch", "dtw_ua_c: Dimensions of input s and t must match.");
    }
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1, nt, mxREAL);
    
    /*  create a C pointer to a copy of the output matrix */
    dp = mxGetPr(plhs[0]);
    
    /*  call the C subroutine */
    dtw_ua_c(s, t, ns, nt, k, dp);
}
