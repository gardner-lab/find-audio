% this is a demo showing the use of our dynamic time warping package 
% we provide both Matlab version and C/MEX version
% the C/MEX version is much faster and highly recommended

% clear;clc;close all;

compile;

a = rand(10, 100);
b = rand(10, 500);

tic;
[d1, s1]=dtw_ua_cos(a, b, 1.5);
t1=toc;

tic;
[d2, s2]=dtw_ua_cos_c(a, b, 1.5);
t2=toc;

fprintf('Using Matlab version: distance=%f, running time=%f\n',min(d1),t1);
fprintf('Using C/MEX version: distance=%f, running time=%f\n',min(d2),t2);

% checks
disp(any(s1(:) ~= s2(:)));
disp(max(abs(d1 - d2)));
ax1 = subplot(2, 1, 1); plot(1:length(d1), d1, 1:length(d1), abs(size(a, 2) - ((1:length(d1)) - s1)));
ax2 = subplot(2, 1, 2); plot(1:length(d2), d2, 1:length(d2), abs(size(a, 2) - ((1:length(d2)) - s2)));
linkaxes([ax1, ax2], 'x');

% benchmark
tms = zeros(1, 1000);
for i = 1:1000
    tic;
    d2=dtw_ua_c(a,b);
    t2=toc;
    tms(i) = t2;
end
disp(mean(tms));
