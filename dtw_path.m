function path = dtw_path(s, t, alphas)
% s: signal 1, size is ns*k, column for time, row for channel 
% t: signal 2, size is nt*k, column for time, row for channel 
% alpha: non-diagonal penalty multiplier
% d: resulting distances

ns = size(s, 2);
nt = size(t, 2);
if size(s, 1) ~= size(t, 1)
    error('Error in dtw the dimensions of the two input signals do not match.');
end
if ~exist('alphas', 'var') || isempty(alphas)
    alphas = 1;
end
if isscalar(alphas)
    alphas = alphas * ones(1, ns);
end
if length(alphas) ~= ns
    error('The length of alpha does not match input signal s.');
end

%% initialization
D = zeros(ns + 1, nt + 1) + Inf; % cache matrix
D(1, 1) = 0; % unanchored

% track steps
steps = zeros(ns + 1, nt + 1);

%% begin dynamic programming
for i = 1:ns
    for j = 1:nt
        % current cost
        oost = norm(s(:, i) - t(:, j));
        
        % transition cost and index
        [cost, idx] = min([D(i, j) + oost, ...
            D(i + 1, j) + oost * alphas(i), ...
            D(i, j + 1) + oost * alphas(i)]);
        
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
