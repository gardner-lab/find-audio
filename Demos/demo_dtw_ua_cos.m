% Demonstrates DTW_UA_COS using two random signals.

% must change diretory to have access to private methods
old_dir = pwd;
new_dir = fileparts(mfilename('fullpath'));
cd(fullfile(new_dir, '..', 'private'));

% generate signal
a = rand(10, 100);
b = rand(10, 500);

% run
tic;
[d, s] = dtw_ua_cos(a, b, 1.5);
t = toc;

% print results
fprintf('distance=%f, running time=%f\n',min(d),t);

% checks
plot(1:length(d), d, 1:length(d), abs(size(a, 2) - ((1:length(d)) - s)));

% benchmark
tms = zeros(1, 1000);
for i = 1:length(tms)
    tic;
    d = dtw_ua(a,b);
    t = toc;
    tms(i) = t;
end
fprintf('Mean time: %f\n', mean(tms));

% return to old directory
cd(old_dir);
