function [scores, starts] = dtw_ua_cos(s, t, alphas)
% s: signal 1, size is ns*k, column for time, row for channel 
% t: signal 2, size is nt*k, column for time, row for channel 
% alpha: non-diagonal penalty multiplier
% d: resulting distances

ns = size(s, 2);
nt = size(t, 2);
if size(s, 1) ~= size(t, 1)
    error('Error in dtw_ua: the dimensions of the two input signals do not match.');
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
D(1, :) = 0; % unanchored

starts = zeros(ns + 1, nt + 1);
starts(1, :) = 0:nt;

%% begin dynamic programming
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
