mex CFLAGS='$FLAGS -Weverything' LDFLAGS='$LDFLAGS -framework Accelerate' dtw_ua_c.c;
mex CFLAGS='$FLAGS -Weverything' LDFLAGS='$LDFLAGS -framework Accelerate' dtw_ua_cos_c.c;
mex CFLAGS='$FLAGS -Weverything' dtpa_c.c;
mex CFLAGS='$FLAGS -Weverything' LDFLAGS='$LDFLAGS -framework Accelerate' dtw_path_c.c;

% mex -DCALCULATE_PATH LDFLAGS='$LDFLAGS -framework Accelerate' dtw_ua_xcorr_c.c;