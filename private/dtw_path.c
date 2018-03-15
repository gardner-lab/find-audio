/**
 * This is the C version of the dtw_path.m function, which offers substantial performance 
 * benefits in terms of faster execution and lower memory usage.
 * 
 * To compile, run the command: `compile_dt_mex`.
 * 
 * Usage mirrors the dtw_path.m function. Run `help dtw_path` for more information.
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

enum step {
    DIAGONAL,
    UP,
    LEFT
};

static mxArray *dtw_path_c(double *mat_template, double *mat_signal, size_t cols_template, size_t cols_signal, size_t rows, double *alphas) {
    /* memory */
    double *mat_score;
    enum step *mat_step;
    
    /* iteration variables */
    size_t i_template, i_signal;
    size_t col_cur, col_last;
    size_t row_cur, row_last;
    
    /* per iteration variables */
    double cost;
    double alpha;
    double t_path, b_path;
    enum step b_step;
    
    /* path building variables */
    double *mat_path;
    size_t step_count;
    
    /* return building variables */
	mxArray *ret;
    double *ptr;
    
    /* allocate memory */
    mat_score = (double *)mxCalloc((cols_template + 1) * (cols_signal + 1), sizeof(double));
    mat_step = (enum step *)mxCalloc((cols_template + 1) * (cols_signal + 1), sizeof(enum step));
    
    /* seed memory */
    for (i_signal = 1; i_signal <= cols_signal; ++i_signal) {
        mat_score[i_signal * (cols_template + 1)] = DBL_MAX;
    }
    for (i_template = 1; i_template <= cols_template; ++i_template) {
        mat_score[i_template] = DBL_MAX;
    }
    
    /* dynamic programming */
    /* for each column of the signal... */
    for (i_signal = 1; i_signal <= cols_signal; ++i_signal) {
        /* iteration variables (easier lookup in matrix) */
        col_cur = i_signal * (cols_template + 1);
        col_last = (i_signal - 1) * (cols_template + 1);
        
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
                b_path = mat_score[col_last + row_last];
                b_step = DIAGONAL;
            }
            else {
                /* diagonal */
                b_path = mat_score[col_last + row_last] + cost;
                b_step = DIAGONAL;
                
                /* up */
                t_path = mat_score[col_cur + row_last] + cost * alpha;
                if (t_path < b_path) {
                    b_path = t_path;
                    b_step = UP;
                }
                
                /* left */
                t_path = mat_score[col_last + row_cur] + cost * alpha;
                if (t_path < b_path) {
                    b_path = t_path;
                    b_step = LEFT;
                }
            }
            
            /* store values */
            mat_score[col_cur + row_cur] = b_path;
            mat_step[col_cur + row_cur] = b_step;
        }
    }
    
    /* add up steps */
    /* reuse mat_score meory to hold path (since doubles are easier to pass back to matlab) */
    mat_path = mat_score + (cols_template + 1) * (cols_signal + 1);
    step_count = 0;
    i_signal = cols_signal;
    i_template = cols_template;
    while (i_template > 0 || i_signal > 0) {
        /* write it */
        mat_path -= 2;
        mat_path[0] = (double)i_template;
        mat_path[1] = (double)i_signal;
        ++step_count;
        
        /* step backwards */
        col_cur = i_signal * (cols_template + 1);
        row_cur = i_template;
        switch (mat_step[col_cur + row_cur]) {
            case DIAGONAL:
                i_signal = i_signal - 1;
                i_template = i_template - 1;
                break;
            case UP:
                i_template = i_template - 1;
                break;
            case LEFT:
                i_signal = i_signal - 1;
                break;
        }
    }
    
    /* generate memory to return */
    ret = mxCreateDoubleMatrix(2, step_count, mxREAL);
    ptr = mxGetPr(ret);
    memcpy(ptr, mat_path, 2 * step_count * sizeof(double));
    
    /* free memory */
    mxFree(mat_score);
    mxFree(mat_step);
    
    return ret;
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
    
    /*  check for proper number of arguments */
    if (nrhs != 2 && nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:dtw_path_c:invalidNumInputs", "Two inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt( "MATLAB:dtw_path_c:invalidNumOutputs", "One output required.");
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
        mexErrMsgIdAndTxt("MATLAB:dtw_path_c:dimNotMatch", "Dimensions of input s and t must match.");
    }
    
    /* get alpha */
    if (nrhs >= 3) {
        /* accept vector or scalar */
        if (mxGetN(prhs[2]) == ns && mxGetM(prhs[2]) == 1) {
            alphas = mxGetPr(prhs[2]);
        }
        else {
            alpha = getScalar(prhs[2], "MATLAB:dtw_path_c:alphaNotScalar", "Alpha must be a scalar or a vector with length matching input s.");
            
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
    
    /*  call the C subroutine */
    plhs[0] = dtw_path_c(s, t, ns, nt, k, alphas);
    
    /* release any memory allocated */
    if (free_alphas) {
        mxFree(alphas);
    }
}
