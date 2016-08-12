function path = dtw_path(s, t, alphas)
%DTW_PATH Dynamic time warping returning best path between two signals
%   Uses traditional dynamic time warping. Scoring is performed using the 
%   L2 norm, and an optional alpha parameter provides either a constant or 
%   a time-varying non-diagonal penalty. The method keeps track of what
%   transition is made at each step of the dynamic programming process, and
%   these steps are retraced at the end to figure out the optimal path.
%
%   Note: There is a C MEX version of this function that is about 100x
%   faster and should be used instead. If you see a warning, run 
%   "compile_dt_mex" once after installation to compile MEX versions.
%
%   PATH = dtw_path(S, T) runs traditional dynamic time warping over
%   signals S and T. S and T must have the same number of rows,
%   representing the features being compared. Columns represent time. PATH
%   will be a 2xN matrix containing corresponding time indices from both S
%   and T resulting from the warping process.
%
%   PATH = dtw_path(S, T, ALPHAS) adds a non-diagonal penalty (as a 
%   multiplier) to the scoring process. ALPHAS can either be a scalar (a 
%   constant penalty throughout the template) or can be a vector with the
%   same length as the number of columns in S, allowing the penalty to vary
%   throughout the template.
%
%   Scoring: L2 norm of the difference between a column in TEMPLATE and a
%   column in SIGNAL, multiplied by the corresponding ALPHAS (or 1) if the
%   last step was not diagonal.

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
D(1, 1) = 0; % anchored

% track steps
steps = zeros(ns + 1, nt + 1);

%% dynamic programming
for i = 1:ns
    for j = 1:nt
        % current cost
        oost = norm(s(:, i) - t(:, j));
        
        % transition cost and index
        if isnan(oost)
            cost = D(i, j);
            idx = 1;
        else
            [cost, idx] = min([D(i, j) + oost, ...
                D(i + 1, j) + oost * alphas(i), ...
                D(i, j + 1) + oost * alphas(i)]);
        end
            
        
        % store
        D(i + 1, j + 1) = cost;
        steps(i + 1, j + 1) = idx;
    end
end

%% assemble path
path = zeros(2, ns + nt + 1);
i = ns;
j = nt;
k = 1;
while i > 0 || j > 0
    path(:, k) = [i; j];
    switch steps(i + 1, j + 1)
        case 1
            i = i - 1;
            j = j - 1;
        case 2
            j = j - 1;
        case 3
            i = i - 1;
    end
    k = k + 1;
end
path = path(:, (k - 1):-1:1);

end
