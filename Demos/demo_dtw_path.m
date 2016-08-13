% Demonstrates DTW_PATH using two random signals.

% must change diretory to have access to private methods
old_dir = pwd;
new_dir = fileparts(mfilename('fullpath'));
cd(fullfile(new_dir, '..', 'private'));

% generate signal
a = rand(10, 100);
b = rand(10, 500);

% run
tic;
p = dtw_path(a, b, 1.5);
t = toc;

% print results
fprintf('steps=%d, running time=%f\n',size(p, 2),t);

% benchmark
tms = zeros(1, 1000);
for i = 1:length(tms)
    tic;
    p = dtw_path(a,b);
    t = toc;
    tms(i) = t;
end
fprintf('Mean time: %f\n', mean(tms));

% return to old directory
cd(old_dir);
