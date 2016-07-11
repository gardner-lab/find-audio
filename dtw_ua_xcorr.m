function [d, p] = dtw_ua_xcorr(s, t, max_lag, alpha)
% s: signal 1, size is ns*k, column for time, row for channel 
% t: signal 2, size is nt*k, column for time, row for channel 
% d: resulting distances

ns = size(s, 2);
nt = size(t, 2);
if size(s, 1) ~= size(t, 1)
    error('Error in dtw_ua: the dimensions of the two input signals do not match.');
end
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 1;
end

% default max lag
if ~exist('max_lag', 'var') || isempty(max_lag)
    max_lag = size(s, 1) - 1;
end

%% initialization
D = zeros(ns + 1, nt + 1) + Inf; % cache matrix
D(1, :) = 0; % unanchored

Pup = zeros(ns + 1, nt + 1, 'int32');
Pleft = zeros(ns + 1, nt + 1, 'int32');

%% begin dynamic programming
for i = 1:ns
    for j = 1:nt
        xc = xcorr(s(:, i), t(:, j), max_lag);
        oost = 1 / max(xc);
        [cost, idx] = min([D(i, j) + oost, ...
            D(i + 1, j) + oost * alpha, ...
            D(i, j + 1) + oost * alpha]);
        
        D(i + 1, j + 1) = cost;
        switch idx
            case 1 % diagonal
                Pleft(i + 1, j + 1) = Pleft(i, j);
                Pup(i + 1, j + 1) = Pup(i, j);
            case 2 % left
                Pleft(i + 1, j + 1) = Pleft(i + 1, j) + 1;
                Pup(i + 1, j + 1) = Pup(i + 1, j);
            case 3 % up
                Pleft(i + 1, j + 1) = Pleft(i, j + 1);
                Pup(i + 1, j + 1) = Pup(i, j + 1) + 1;
        end
    end
end

%% output
d = D(end, 2:end);
p = [Pup(end, 2:end); Pleft(end, 2:end)];

end
