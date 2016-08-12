function [scores, start] = dtpa(s, t, max_lag, alphas)
%DTPA Unanchored dynamic time warping with translation
%   Uses an unanchored version of dynamic time warping to search for a
%   template in a signal, allowing translations across the feature vectors.
%   Unlike traditional dynamic time warping, this does not anchor the 
%   beginning and end of the two data streams, but rather looks for any 
%   arbitrarily occurence of the template in the signal. In addition, the
%   dynamic programming is constructed as a 3D problem where each step will
%   move either through the signal or the template, as well as through 
%   different lags (or translations). This is useful, for example, if
%   matching a log scaled spectrogram and there may be minor pitch
%   variations. The amount of translation can only increase or decrease by
%   1 step at each time step. Scoring is performed using the L2 norm, and 
%   an optional alpha parameter provides either a constant or a 
%   time-varying non-diagonal penalty.
%
%   Note: There is a C MEX version of this function that is about 100x
%   faster and should be used instead. If you see a warning, run 
%   "compile_dt_mex" once after installation to compile MEX versions.
%
%   SCORES = dtpa(TEMPLATE, SIGNAL) performs unanchored dynamic time look
%   for an arbitrarily positioned occurence of TEMPLATE in SIGNAL. TEMPLATE
%   and SIGNAL must have the same number of rows, representing the features
%   being compared. Columns represent time. SCORES will return a vector of
%   scores with the same length as the number of columns in SIGNAL
%   representing the best match ending at each step of the signal.
%
%   SCORES = dtpa(TEMPLATE, SIGNAL, MAX_LAG) sets the amount of translation
%   that can occur between the feature vectors. If ommitted, the default
%   values allows offsets up to half the number of features.
%
%   SCORES = dtpa(TEMPLATE, SIGNAL, MAX_LAG, ALPHAS) adds a non-diagonal 
%   penalty (as a multiplier) to the scoring process. ALPHAS can either be 
%   a scalar (a constant penalty throughout the template) or can be a 
%   vector with the same length as the number of columns in TEMPLATE,
%   allowing the penalty to vary throughout the template.
%
%   [SCORES, STARTS] = dtpa(...) returns the STARTS time step that
%   corresponds to the end of each match. This allows figuring out the
%   start and end point corresponding to each score.
%
%   Scoring: L2 norm of the difference between a column in TEMPLATE and a
%   column in SIGNAL, where values are offset by the lag parameter.
%   After translation, features without a corresponding match are ignored.
%   The score at each step is multiplied by the corresponding ALPHAS (or 1)
%   if the last step was not diagonal.

% display warning
warning('Non-MEX dynamic time warping being used. Call "compile_dt_mex" once to generate faster MEX functions.');

%% validate inputs

% sizes
nv = size(s, 1);
ns = size(s, 2);
nt = size(t, 2);
if ~exist('max_lag', 'var') || isempty(max_lag)
    max_lag = round(size(s, 1) / 2);
end
if nv ~= size(t, 1)
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

%% pad
if 0 < max_lag
    s = [zeros(max_lag, ns); s; zeros(max_lag, ns)];
end

%% setup
% maximum offset
max_offset = 1 + 2 * max_lag;

% preallocate matrix
D = zeros(ns + 1, nt + 1, max_offset) + Inf;
D(1, :, :) = 0; % unanchored

start = zeros(ns + 1, nt + 1, 1 + 2 * max_lag);
start(1, :, :) = repmat(0:nt, 1, 1, 1 + 2 * max_lag);

%% dynamic programming
for i = 1:ns
    for j = 1:nt
        for k = 1:max_offset
            % shift norm difference
            oost = norm(s(k:(k + nv - 1), i) - t(:, j));
            
            % indices for potential translation
            if k > 1 && k < max_offset
                r = [k k - 1 k + 1];
            elseif k > 1
                r = [k k - 1];
            else
                r = [k k + 1];
            end
            
            % transitions
            transitions = [...
                squeeze(D(i, j, r)) + oost ... % diagonals
                squeeze(D(i + 1, j, r)) + oost * alphas(i) ... % left
                squeeze(D(i, j + 1, r)) + oost * alphas(i)... % up
                ]';
            
            starts = [
                squeeze(start(i, j, r)) ... % diagonals
                squeeze(start(i + 1, j, r)) ... % left
                squeeze(start(i, j + 1, r))... % up
                ]';
            
            
            % find smallest transition cost
            [cost, idx] = min(transitions(:));
            D(i + 1, j + 1, k) = cost;
            start(i + 1, j + 1, k) = starts(idx);
        end
    end
end

%% output
[scores, idx] = min(D(end, 2:end, :), [], 3);
start = start(sub2ind(size(start), ns + ones(1, 500), 2:(nt + 1), idx));

end
