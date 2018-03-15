/**
 * This is the C version of the dtw_ua.m function, which offers substantial performance 
 * benefits in terms of faster execution and lower memory usage.
 * 
 * To compile, run the command: `compile_dt_mex`.
 * 
 * Usage mirrors the dtw_ua.m function. Run `help dtw_ua` for more information.
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

static double vectorDistance(double *s, double *t, size_t k) {
    double result = 0;
    
#if __APPLE__
    /* surprisingly, this may not be faster... */
    vDSP_distancesqD(s, 1, t, 1, &result, k);
#else
    double ss, tt;
    size_t x;
    for (x = 0; x < k; x++) {
        ss = s[x];
        tt = t[x];
        result += ((ss - tt) * (ss - tt));
    }
#endif
    
    result = sqrt(result);
    
    return result;
}

static void dtw_ua_c(double *mat_template, double *mat_signal, size_t cols_template, size_t cols_signal, size_t rows, double *alphas, double *out_score, double *out_start) {
    /* memory */
    double *mat_score[2]; /* matrix of scores */
    size_t *mat_start[2]; /* matrix of starts */
    
    /* iteration variables */
    size_t i_template, i_signal;
    size_t col_cur, col_last;
    size_t row_cur, row_last;
    
    /* per iteration variables */
    double cost;
    double alpha;
    double t_path, b_path;
    size_t b_start;
    
    /* allocate memory */
    mat_score[0] = (double *)mxCalloc(cols_template + 1, sizeof(double));
    mat_score[1] = (double *)mxCalloc(cols_template + 1, sizeof(double));
    
    mat_start[0] = (size_t *)mxCalloc(cols_template + 1, sizeof(size_t));
    mat_start[1] = (size_t *)mxCalloc(cols_template + 1, sizeof(size_t));
    
    /* seed memory */
    for (i_template = 1; i_template <= cols_template; ++i_template) {
        row_cur = i_template;
        mat_score[0][row_cur] = DBL_MAX;
    }
    
    /* dynamic programming */
    /* for each column of the signal... */
    for (i_signal = 1; i_signal <= cols_signal; ++i_signal) {
        /* iteration variables (easier lookup in matrix) */
        col_cur = i_signal % 2;
        col_last = (i_signal - 1) % 2;
        
        /* seed start */
        mat_start[col_cur][0] = i_signal;
        
        /* for each column of the tempalte... */
        for (i_template = 1; i_template <= cols_template; ++i_template) {
            /* iteration variables (easier lookup in matrix) */
            row_cur = i_template;
            row_last = i_template - 1;
            
            /* get current alpha */
            alpha = alphas[i_template - 1];
            
            /* calculate norm */
            cost = vectorDistance(mat_signal + rows * (i_signal - 1), mat_template + rows * (i_template - 1), rows);
            
            /* special is nan */
            if isnan(cost) {
                /* assume diagonal */
                b_path = mat_score[col_last][row_last];
                b_start = mat_start[col_last][row_last];
            }
            else {
                /* diagonal */
                b_path = mat_score[col_last][row_last] + cost;
                b_start = mat_start[col_last][row_last];
                
                /* up */
                t_path = mat_score[col_cur][row_last] + cost * alpha;
                if (t_path < b_path) {
                    b_path = t_path;
                    b_start = mat_start[col_cur][row_last];
                }
                
                /* left */
                t_path = mat_score[col_last][row_cur] + cost * alpha;
                if (t_path < b_path) {
                    b_path = t_path;
                    b_start = mat_start[col_last][row_cur];
                }
            }
            
            /* store values */
            mat_score[col_cur][row_cur] = b_path;
            mat_start[col_cur][row_cur] = b_start;
        }
        
        /* store values */
        row_cur = cols_template; /* does not actually need to be reset */
        out_score[i_signal - 1] = mat_score[col_cur][row_cur];
        if (out_start) {
            out_start[i_signal - 1] = (double)mat_start[col_cur][row_cur];
        }
    }
    
    /* free memory */
    mxFree(mat_score[0]);
    mxFree(mat_score[1]);
    mxFree(mat_start[0]);
    mxFree(mat_start[1]);
}

static double getScalar(const mxArray *in, const char *err_id, const char *err_str) {
    /* check scalar */
    if (!mxIsDouble(in) || mxIsComplex(in) || mxGetN(in) * mxGetM(in) != 1) {
        mexErrMsgIdAndTxt(err_id, err_str);
    }
    
    /* get the scalar input */
    return mxGetScalar(in);
}

/* the gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *s, *t;
    size_t ns, nt, k;
    size_t i;
    double alpha;
    double *alphas;
    bool free_alphas = false;
    double *dp;
    double *dq;
    
    /*  check for proper number of arguments */
    if (nrhs != 2 && nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:dtw_ua_c:invalidNumInputs", "Two inputs required.");
    }
    if (nlhs != 1 && nlhs != 2) {
        mexErrMsgIdAndTxt( "MATLAB:dtw_ua_c:invalidNumOutputs", "One or two outputs required.");
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
        mexErrMsgIdAndTxt("MATLAB:dtw_ua_c:dimNotMatch", "Dimensions of input s and t must match.");
    }
    
    /* get alpha */
    if (nrhs >= 3) {
        /* accept vector or scalar */
        if (mxGetN(prhs[2]) == ns && mxGetM(prhs[2]) == 1) {
            alphas = mxGetPr(prhs[2]);
        }
        else {
            alpha = getScalar(prhs[2], "MATLAB:dtw_ua_c:alphaNotScalar", "Alpha must be a scalar or a vector with length matching input s.");
            
            free_alphas = true;
            alphas = (double *)mxMalloc(ns * sizeof(double));
            for (i = 0; i < ns; ++i) {
                alphas[i] = alpha;
            }
        }
    }
    else {
        free_alphas = true;
        alphas = (double *)mxMalloc(ns * sizeof(double));
        for (i = 0; i < ns; ++i) {
            alphas[i] = 1;
        }
    }
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1, nt, mxREAL);
    
    /*  create a C pointer to a copy of the output matrix */
    dp = mxGetPr(plhs[0]);
    
    /* save path information */
    if (nlhs == 2) {
        /* create output matrix */
        plhs[1] = mxCreateDoubleMatrix(1, nt, mxREAL);
        
        /* create a C pointer */
        dq = mxGetPr(plhs[1]);
    }
    else {
        dq = NULL;
    }
    
    /*  call the C subroutine */
    dtw_ua_c(s, t, ns, nt, k, alphas, dp, dq);
    
    /* release any memory allocated */
    if (free_alphas) {
        mxFree(alphas);
    }
}
