% Copyright (C) 2013 Quan Wang <wangq10@rpi.edu>,
% Signal Analysis and Machine Perception Laboratory,
% Department of Electrical, Computer, and Systems Engineering,
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

% this is a demo showing the use of our dynamic time warping package 
% we provide both Matlab version and C/MEX version
% the C/MEX version is much faster and highly recommended

% clear;clc;close all;

compile;

a=rand(10, 100);
b=rand(10, 500);

tic;
d1=dtw_ua(a,b);
t1=toc;

tic;
d2=dtw_ua_c(a,b);
t2=toc;

fprintf('Using Matlab version: distance=%f, running time=%f\n',min(d1),t1);
fprintf('Using C/MEX version: distance=%f, running time=%f\n',min(d2),t2);

% checks
disp(mean(d1 - d2));
subplot(2,1,1);plot(d1);
subplot(2,1,2);plot(d2);


tms = zeros(1, 1000);
for i = 1:1000
    tic;
    d2=dtw_ua_c(a,b);
    t2=toc;
    tms(i) = t2;
end
disp(mean(tms));

