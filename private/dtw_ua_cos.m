function [scores, starts] = dtw_ua_cos(s, t, alphas)
%DTW_UA_COS Unanchored dynamic time warping to find template in signal
%   Uses an unanchored version of dynamic time warping to search for a
%   template in a signal. Unlike traditional dynamic time warping, this
%   does not anchor the beginning and end of the two data streams, but
%   rather looks for any arbitrarily occurence of the template in the
%   signal. Scoring is performed using the cosine similarity, and an 
%   optional alpha parameter provides either a constant or a time-varying
%   non-diagonal penalty.
%
%   Note: There is a C MEX version of this function that is about 100x
%   faster and should be used instead. If you see a warning, run 
%   "compile_dt_mex" once after installation to compile MEX versions.
%
%   SCORES = dtw_ua_cos(TEMPLATE, SIGNAL) performs unanchored dynamic time 
%   look for an arbitrarily positioned occurence of TEMPLATE in SIGNAL. 
%   TEMPLATE and SIGNAL must have the same number of rows, representing the 
%   features being compared. Columns represent time. SCORES will return a 
%   vector of scores with the same length as the number of columns in 
%   SIGNAL representing the best match ending at each step of the signal.
%
%   SCORES = dtw_ua_cos(TEMPLATE, SIGNAL, ALPHAS) adds a non-diagonal 
%   penalty (as a multiplier) to the scoring process. ALPHAS can either be 
%   a scalar (a constant penalty throughout the template) or can be a 
%   vector with the same length as the number of columns in TEMPLATE, 
%   allowing the penalty to vary throughout the template.
%
%   [SCORES, STARTS] = dtw_ua_cos(...) returns the STARTS time step that
%   corresponds to the end of each match. This allows figuring out the
%   start and end point corresponding to each score.
%
%   Scoring: The cosine between a column in TEMPLATE and a column in 
%   SIGNAL, multiplied by the corresponding ALPHAS (or 1) if the last step
%   was not diagonal.

% display warning
warning('Non-MEX dynamic time warping being used. Call "compile_dt_mex" once to generate faster MEX functions.');

%% validate inputs

% sizes
ns = size(s, 2);
nt = size(t, 2);
if size(s, 1) ~= size(t, 1)
    error('The number of features in the two input signals do not match.');
end

% alphas
if ~exist('alphas', 'var') || isempty(alphas)
    alphas = 1;
end
if isscalar(alphas)
    alphas = alphas * ones(1, ns);
end
if length(alphas) ~= ns
    error('The length of alpha does not match input template.');
end

%% setup
% preallocate matrix
D = zeros(ns + 1, nt + 1) + Inf;
D(1, :) = 0; % unanchored

starts = zeros(ns + 1, nt + 1);
starts(1, :) = 0:nt;

%% dynamic programming
for i = 1:ns
    for j = 1:nt
        oost = dot(s(:, i), t(:, j)) / (norm(s(:, i)) * norm(t(:, j)));
        [cost, idx] = min([D(i, j) + oost, ...
            D(i + 1, j) + oost * alphas(i), ...
            D(i, j + 1) + oost * alphas(i)]);
        
        D(i + 1, j + 1) = cost;
        switch idx
            case 1 % diagonal
                starts(i + 1, j + 1) = starts(i, j);
            case 2 % left
                starts(i + 1, j + 1) = starts(i + 1, j);
            case 3 % up
                starts(i + 1, j + 1) = starts(i, j + 1);
        end
    end
end

%% output
scores = D(end, 2:end);
starts = starts(end, 2:end);

end
