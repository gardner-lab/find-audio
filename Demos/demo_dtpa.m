% Demonstrates DTPA using two random signals with a section that has a
% the same values but at a slight offset.

% must change diretory to have access to private methods
old_dir = pwd;
new_dir = fileparts(mfilename('fullpath'));
cd(fullfile(new_dir, '..', 'private'));

% generate signal
a = rand(10, 100);
b = rand(10, 500);

% make correlated
b(3:10, 301:400) = a(1:8, :);
a(1:8, :) = a(1:8, :);

% alpha
alphas = 10 * 0.9 .^ ([0:49 49:-1:0]) + 2;

% run
tic;
[d, s] = dtpa(a, b, 3, alphas);
t = toc;

% print results
fprintf('distance=%f, running time=%f\n', min(d), t);

% plot
plot(1:length(d), d, 1:length(d), abs(size(a, 2) - ((1:length(d)) - s)));

% return to directory
cd(old_dir);
