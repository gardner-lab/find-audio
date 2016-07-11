/**
 * This is the C/MEX code of dynamic time warping of two signals, using xcorr as a score
 *
 * compile: 
 *     mex dtw_ua_xcorr_c.c
 *
 * usage:
 *     d=dtw_ua_xcorr_c(s,t[,maxLag])
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

double vectorXcorr(double *s, double *t, int k, int maxLag, double *mem_in, double *mem_out) {
    double result = 0;
    
    // copy into memory in
    memcpy(mem_in + maxLag, s, k * sizeof(double));
    
#if __APPLE__
    vDSP_convD(mem_in, 1, t, 1, mem_out, 1, 1 + 2 * maxLag, k);
    vDSP_maxvD(mem_out, 1, &result, 1 + 2 * maxLag);
    //printf("[%f %f %f ...] xcorr [%f %f %f ...] = %f\n", s[0], s[1], s[2], t[0], t[1], t[2], result);
#else
#error "Write me!"
#endif
    return result;
}

#ifdef CALCULATE_PATH
struct path {
    int up;
    int left;
};
#endif

void dtw_ua_xcorr_c(double *s, double *t, int ns, int nt, int k, int maxLag, double alpha, double *dp, int *dq) {
    double *xcorr_mem_in;
    double *xcorr_mem_out;
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
    
    // xcorr setup
    if (maxLag == -1 || maxLag > (k - 1)) {
        maxLag = k - 1;
    }
    xcorr_mem_in = (double *)mxCalloc(k + (maxLag * 2), sizeof(double));
    xcorr_mem_out = (double *)mxCalloc(1 + (maxLag * 2), sizeof(double));
    
    // dynamic programming
    for (i = 1; i <= ns; i++) {
        cur = i % 2;
        last = (i - 1) % 2;
        
        D[cur][0] = DBL_MAX;
        
        for (j = 1; j <= nt; j++) {
            // calculate norm
            cost = 1 / vectorXcorr(s + k * (i - 1), t + k * (j - 1), k, maxLag, xcorr_mem_in, xcorr_mem_out);
            
            a = D[last][j] + cost * alpha; // up
            b = D[cur][j - 1] + cost * alpha; // left
            c = D[last][j - 1] + cost; // diagonal
            
            if (a < b && a < c) {
                // up
                D[cur][j] = a;
#ifdef CALCULATE_PATH
                P[cur][j] = P[last][j];
                P[cur][j].up++;
#endif
            }
            else if (b < c) {
                // left
                D[cur][j] = b;
#ifdef CALCULATE_PATH
                P[cur][j] = P[cur][j - 1];
                P[cur][j].left++;
#endif
            }
            else {
                // diagonal
                D[cur][j] = c;
#ifdef CALCULATE_PATH
                P[cur][j] = P[last][j - 1];
#endif
            }
        }
    }
    
    // xcorr cleanup
    mxFree(xcorr_mem_in);
    mxFree(xcorr_mem_out);
    
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

double getScalar(const mxArray *in, const char *err_id, const char *err_str) {
    /* check scalar */
    if (!mxIsDouble(in) || mxIsComplex(in) || mxGetN(in) * mxGetM(in) != 1) {
        mexErrMsgIdAndTxt(err_id, err_str);
    }
    
    /* get the scalar input */
    return mxGetScalar(in);
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *s,*t;
    int w;
    int ns, nt, k;
    int maxLag;
    double alpha;
    double *dp;
    int *dq;
    
    /*  check for proper number of arguments */
    if (nrhs < 2 || nrhs > 4) {
        mexErrMsgIdAndTxt("MATLAB:dtw_ua_xcorr_c:invalidNumInputs", "Two to four inputs required.");
    }
#ifdef CALCULATE_PATH
    if (nlhs != 1 && nlhs != 2) {
        mexErrMsgIdAndTxt( "MATLAB:dtw_ua_xcorr_c:invalidNumOutputs", "One or two outputs required.");
    }
#else
    if (nlhs != 1) {
        mexErrMsgIdAndTxt( "MATLAB:dtw_ua_xcorr_c:invalidNumOutputs", "One output required.");
    }
#endif
    
    /* get maxLag */
    if (nrhs >= 3) {
        maxLag = (int)getScalar(prhs[2], "MATLAB:dtw_ua_xcorr_c:maxLagNotScalar", "Max lag must be a scalar.");
    }
    else {
        maxLag = -1;
    }
    
    /* get alpha */
    if (nrhs >= 4) {
        alpha = getScalar(prhs[3], "MATLAB:dtw_ua_xcorr_c:alphaNotScalar", "Alpha must be a scalar.");
    }
    else {
        alpha = 1;
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
        mexErrMsgIdAndTxt("MATLAB:dtw_ua_xcorr_c:dimNotMatch", "Dimensions of input s and t must match.");
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
    dtw_ua_xcorr_c(s, t, ns, nt, k, maxLag, alpha, dp, dq);
}
