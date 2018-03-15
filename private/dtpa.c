/**
 * This is the C version of the dtpa.m function, which offers substantial performance 
 * benefits in terms of faster execution and lower memory usage.
 * 
 * To compile, run the command: `compile_dt_mex`.
 * 
 * Usage mirrors the dtpa.m function. Run `help dtpa` for more information.
 */

#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

static double vectorDistance(double *s, double *t, int len, int lag) {
    double result = 0;
    double ss, tt; /* hold current value from vector */
    int x, y; /* indices into vector */
    
    for (x = 0; x < len; x++) {
        /* apply lag to index */
        y = x + lag;
        
        /* values */
        ss = s[x];
        tt = (y >= 0 && y < len ? t[y] : 0);
        
        result += ((ss - tt) * (ss - tt));
    }
    
    result = sqrt(result);
    
    return result;
}

static void dtw_v2_c(double *mat_template, double *mat_signal, size_t cols_template, size_t cols_signal, size_t rows, size_t max_lag, double *alphas, double *out_score, double *out_start) {
    /* memory */
    double *mat_score[2]; /* matrix of scores */
    size_t *mat_start[2]; /* matrix of starts */
    
    /* iteration variables */
    size_t i_template, i_signal, i_offset;
    size_t num_offsets;
    size_t col_cur, col_last;
    size_t row_cur, row_last;
    
    /* per iteration variables */
    double cost;
    double alpha;
    double t_path, b_path;
    size_t b_start;
    
    /* lag setup */
    num_offsets = 1 + 2 * max_lag;
    
    /* allocate memory */
    mat_score[0] = (double *)mxCalloc((cols_template + 1) * num_offsets, sizeof(double));
    mat_score[1] = (double *)mxCalloc((cols_template + 1) * num_offsets, sizeof(double));
    
    mat_start[0] = (size_t *)mxCalloc((cols_template + 1) * num_offsets, sizeof(size_t));
    mat_start[1] = (size_t *)mxCalloc((cols_template + 1) * num_offsets, sizeof(size_t));
    
    /* seed memory */
    for (i_template = 1; i_template <= cols_template; ++i_template) {
        row_cur = i_template * num_offsets;
        for (i_offset = 0; i_offset < num_offsets; ++i_offset) {
            mat_score[0][row_cur + i_offset] = DBL_MAX;
        }
    }
    
    /* dynamic programming */
    /* for each column of the signal... */
    for (i_signal = 1; i_signal <= cols_signal; ++i_signal) {
        /* iteration variables (easier lookup in matrix) */
        col_cur = i_signal % 2;
        col_last = (i_signal - 1) % 2;
        
        /* seed start */
        for (i_offset = 0; i_offset < num_offsets; ++i_offset) {
            mat_start[col_cur][i_offset] = i_signal;
        }
        
        /* for each column of the tempalte... */
        for (i_template = 1; i_template <= cols_template; ++i_template) {
            /* iteration variables (easier lookup in matrix) */
            row_cur = i_template * num_offsets;
            row_last = (i_template - 1) * num_offsets;
            
            /* get current alpha */
            alpha = alphas[i_template - 1];
            
            /* for each offset... */
            for (i_offset = 0; i_offset < num_offsets; ++i_offset) {
                /* calculate cost */
                cost = vectorDistance(mat_signal + rows * (i_signal - 1), mat_template + rows * (i_template - 1), (int)rows, (int)i_offset - (int)max_lag);
                
                /* diagonal, same offset */
                b_path = mat_score[col_last][row_last + i_offset] + cost;
                b_start = mat_start[col_last][row_last + i_offset];
                
                /* up, same offset */
                t_path = mat_score[col_cur][row_last + i_offset] + cost * alpha;
                if (t_path < b_path) {
                    b_path = t_path;
                    b_start = mat_start[col_cur][row_last + i_offset];
                }
                
                /* left, same offset */
                t_path = mat_score[col_last][row_cur + i_offset] + cost * alpha;
                if (t_path < b_path) {
                    b_path = t_path;
                    b_start = mat_start[col_last][row_cur + i_offset];
                }
                
                /* shift offset down */
                if (i_offset > 0) {
                    /* diagonal, shift offset */
                    t_path = mat_score[col_last][row_last + i_offset - 1] + cost;
                    if (t_path < b_path) {
                        b_path = t_path;
                        b_start = mat_start[col_last][row_last + i_offset - 1];
                    }
                    
                    /* up, shift offset */
                    t_path = mat_score[col_cur][row_last + i_offset - 1] + cost * alpha;
                    if (t_path < b_path) {
                        b_path = t_path;
                        b_start = mat_start[col_cur][row_last + i_offset - 1];
                    }
                    
                    /* left, shift offset */
                    t_path = mat_score[col_last][row_cur + i_offset - 1] + cost * alpha;
                    if (t_path < b_path) {
                        b_path = t_path;
                        b_start = mat_start[col_last][row_cur + i_offset - 1];
                    }
                }
                
                /* shift offset up */
                if ((i_offset + 1) < num_offsets) {
                    /* diagonal, shift offset */
                    t_path = mat_score[col_last][row_last + i_offset + 1] + cost;
                    if (t_path < b_path) {
                        b_path = t_path;
                        b_start = mat_start[col_last][row_last + i_offset + 1];
                    }
                    
                    /* up, shift offset */
                    t_path = mat_score[col_cur][row_last + i_offset + 1] + cost * alpha;
                    if (t_path < b_path) {
                        b_path = t_path;
                        b_start = mat_start[col_cur][row_last + i_offset + 1];
                    }
                    
                    /* left, shift offset */
                    t_path = mat_score[col_last][row_cur + i_offset + 1] + cost * alpha;
                    if (t_path < b_path) {
                        b_path = t_path;
                        b_start = mat_start[col_last][row_cur + i_offset + 1];
                    }
                }
                
                /* store values */
                mat_score[col_cur][row_cur + i_offset] = b_path;
                mat_start[col_cur][row_cur + i_offset] = b_start;
            }
        }
        
        /* store values */
        row_cur = cols_template * num_offsets; /* does not actually need to be reset */
        /* POTENTIALLY ACCELERATE */
        for (i_offset = 0; i_offset < num_offsets; ++i_offset) {
            if (0 == i_offset || mat_score[col_cur][row_cur + i_offset] < out_score[i_signal - 1]) {
                out_score[i_signal - 1] = mat_score[col_cur][row_cur + i_offset];
                if (out_start) {
                    out_start[i_signal - 1] = (double)mat_start[col_cur][row_cur + i_offset];
                }
            }
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
    size_t max_lag;
    double alpha;
    double *alphas;
    bool free_alphas = false;
    double *dp;
    double *dq;
    
    /*  check for proper number of arguments */
    if (nrhs < 2 || nrhs > 4) {
        mexErrMsgIdAndTxt("MATLAB:dtw_v2_c:invalidNumInputs", "Two to four inputs required.");
    }
    if (nlhs != 1 && nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:dtw_v2_c:invalidNumOutputs", "One or two outputs required.");
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
        mexErrMsgIdAndTxt("MATLAB:dtw_v2_c:dimNotMatch", "Dimensions of input s and t must match.");
    }
    
    /* get maxLag */
    if (nrhs >= 3) {
        max_lag = (size_t)getScalar(prhs[2], "MATLAB:dtw_v2_c:maxLagNotScalar", "Max lag must be a scalar.");
        
        // bound it by 1 less than the number of rows
        if (max_lag > (k - 1)) {
        	max_lag = k - 1;
        }
    }
    else {
        max_lag = k - 1;
    }
    
    /* get alpha */
    if (nrhs >= 4) {
        /* accept vector or scalar */
        if (mxGetN(prhs[3]) == ns && mxGetM(prhs[3]) == 1) {
            alphas = mxGetPr(prhs[3]);
        }
        else {
            alpha = getScalar(prhs[3], "MATLAB:dtw_v2_c:alphaNotScalar", "Alpha must be a scalar or a vector with length matching input s.");
            
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
    dtw_v2_c(s, t, ns, nt, k, max_lag, alphas, dp, dq);
    
    /* release any memory allocated */
    if (free_alphas) {
        mxFree(alphas);
    }
}
