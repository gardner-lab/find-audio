function [d, start] = dtpa(s, t, max_lag, alphas)
% s: signal 1 (template), size is ns*k, column for time, row for channel 
% t: signal 2 (signal), size is nt*k, column for time, row for channel 
% alphas: non-diagonal penalty multiplier, can be scalar or length equal to
% s
% d: resulting distances

nv = size(s, 1);
ns = size(s, 2);
nt = size(t, 2);
if ~exist('max_lag', 'var') || isempty(max_lag)
    max_lag = round(size(s, 1) / 2);
end
if nv ~= size(t, 1)
    error('The dimensions of the two input signals do not match.');
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

%% pad
if 0 < max_lag
    s = [zeros(max_lag, ns); s; zeros(max_lag, ns)];
end

%% initialization
max_offset = 1 + 2 * max_lag;

D = zeros(ns + 1, nt + 1, max_offset) + Inf; % cache matrix
D(1, :, :) = 0; % unanchored

start = zeros(ns + 1, nt + 1, 1 + 2 * max_lag);
start(1, :, :) = repmat(0:nt, 1, 1, 1 + 2 * max_lag);


%% begin dynamic programming
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
[d, idx] = min(D(end, 2:end, :), [], 3);
start = start(sub2ind(size(start), ns + ones(1, 500), 2:(nt + 1), idx));

end
