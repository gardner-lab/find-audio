function d = dtw_ua(s, t)
% s: signal 1, size is ns*k, column for time, row for channel 
% t: signal 2, size is nt*k, column for time, row for channel 
% d: resulting distances

ns=size(s,2);
nt=size(t,2);
if size(s,1)~=size(t,1)
    error('Error in dtw_ua(): the dimensions of the two input signals do not match.');
end

%% initialization
D=zeros(ns + 1, nt + 1) + Inf; % cache matrix
D(1, :) = 0; % unanchored

%% begin dynamic programming
for i = 1:ns
    for j = 1:nt
        oost = norm(s(:, i) - t(:, j));
        D(i + 1, j + 1) = oost + min([D(i, j + 1) D(i + 1, j) D(i,j)]);
    end
end
d = D(end, 2:end);
