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
#if __APPLE__
#include <Accelerate/Accelerate.h>
#endif

double vectorDistance(double *s, double *t, int k) {
    double result=0;
#if __APPLE__
    // surprisingly, this may not be faster...
    vDSP_distancesqD(s, 1, t, 1, &result, k);
#else
    double ss, tt;
    int x;
    for (x = 0; x < k; x++) {
        ss = s[x];
        tt = t[x];
        result += ((ss - tt) * (ss - tt));
    }
#endif
    result = sqrt(result);
    return result;
}

#ifdef CALCULATE_PATH
struct path {
    int up;
    int left;
};
#endif

void dtw_ua_c(double *s, double *t, int ns, int nt, int k, double *dp, int *dq) {
    double *D[2];
#ifdef CALCULATE_PATH
    struct path *P[2];
#endif
    int i, j;
    int j1, j2;
    double cost;
    double a, b, c;
    int cur, last;
    
    // create D
    D[0] = (double *)mxCalloc(nt + 1, sizeof(double));
    D[1] = (double *)mxCalloc(nt + 1, sizeof(double));
    
#ifdef CALCULATE_PATH
    P[0] = (struct path *)mxCalloc(nt + 1, sizeof(struct path));
    P[1] = (struct path *)mxCalloc(nt + 1, sizeof(struct path));
#endif
    
    // dynamic programming
    for (i = 1; i <= ns; i++) {
        cur = i % 2;
        last = (i - 1) % 2;
        
        D[cur][0] = DBL_MAX;
        
        for (j = 1; j <= nt; j++) {
            // calculate norm
            cost = vectorDistance(s + k * (i - 1), t + k * (j - 1), k);
            
            a = D[last][j]; // up
            b = D[cur][j - 1]; // left
            c = D[last][j - 1]; // diagonal
            
            if (a < b && a < c) {
                // up
                D[cur][j] = cost + a;
#ifdef CALCULATE_PATH
                P[cur][j] = P[last][j];
                P[cur][j].up++;
#endif
            }
            else if (b < c) {
                // left
                D[cur][j] = cost + b;
#ifdef CALCULATE_PATH
                P[cur][j] = P[cur][j - 1];
                P[cur][j].left++;
#endif
            }
            else {
                // diagonal
                D[cur][j] = cost + c;
#ifdef CALCULATE_PATH
                P[cur][j] = P[last][j - 1];
#endif
            }
        }
    }
    
    /* copy to output vector */
    memcpy(dp, D[ns % 2] + 1, nt * sizeof(double));
    
#ifdef CALCULATE_PATH
    /* copy to output vcetor */
    if (dq) {
        memcpy(dq, P[ns % 2] + 1, nt * sizeof(struct path));
    }
#endif
    
    // free memory
    mxFree(D[0]);
    mxFree(D[1]);
#ifdef CALCULATE_PATH
    mxFree(P[0]);
    mxFree(P[1]);
#endif
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *s,*t;
    int w;
    int ns,nt,k;
    double *dp;
    int *dq;
    double *result;
    
    /*  check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:dtw_ua_c:invalidNumInputs", "Two inputs required.");
    }
#ifdef CALCULATE_PATH
    if (nlhs != 1 && nlhs != 2) {
        mexErrMsgIdAndTxt( "MATLAB:dtw_ua_c:invalidNumOutputs", "One or two outputs required.");
    }
#else
    if (nlhs != 1) {
        mexErrMsgIdAndTxt( "MATLAB:dtw_ua_c:invalidNumOutputs", "One output required.");
    }
#endif
    
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
        mexErrMsgIdAndTxt("MATLAB:dtw_ua_c:dimNotMatch", "Dimensions of input s and t must match.");
    }
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1, nt, mxREAL);
    
    /*  create a C pointer to a copy of the output matrix */
    dp = mxGetPr(plhs[0]);
    
#ifdef CALCULATE_PATH
    /* save path information */
    if (nlhs == 2) {
        /* create output matrix */
        plhs[1] = mxCreateNumericMatrix(2, nt, (sizeof(int) == 8 ? mxINT64_CLASS : mxINT32_CLASS), mxREAL);
        
        /* create a C pointer */
        dq = (int *)mxGetPr(plhs[1]);
    }
#endif
    
    /*  call the C subroutine */
    dtw_ua_c(s, t, ns, nt, k, dp, dq);
}
