function spect = normalize_spectrogram(spect, varargin)
%NORMALIZE_SPECTROGRAM Summary of this function goes here
%   Detailed explanation goes here

silence_smooth_win = 5;
silence_smooth_thr = 3;
silence_threshold = [];
spect_log = false;
debug = false;

% load custom parameters
nparams = length(varargin);
if 0 < mod(nparams, 2)
    error('Parameters must be specified as parameter/value pairs');
end
for i = 1:2:nparams
    nm = lower(varargin{i});
    if ~exist(nm, 'var')
        error('Invalid parameter: %s.', nm);
    end
    eval([nm ' = varargin{i+1};']);
end

%% check inputs
if 1 ~= mod(silence_smooth_win, 2)
    error('Silence smoothing must be an odd number.');
end

%% prepare input
if ~isreal(spect)
    spect = abs(spect);
end

%% figure out moments of silence

% figure out threshold
spect_power = sum(log(spect));
if isempty(silence_threshold)
    % two step thresholding: first halfway between max 
    spect_power_med_max = quantile(spect_power, [0.5 0.99]);
    spect_power_thr1 = sum(spect_power_med_max) / 2;
    
    
    spect_power_mn = mean(spect_power(spect_power < spect_power_thr1));
    spect_power_std = std(spect_power(spect_power < spect_power_thr1));
    silence_threshold = spect_power_mn + 3 * spect_power_std;
    
    fprintf('Silence threshold: %f\n', silence_threshold);
    
    % threshold
    if debug
        figure;
        ax1 = subplot(2, 1, 1); imagesc(1:size(spect, 2), 1:size(spect, 1), log(spect));
        ax2 = subplot(2, 1, 2); plot(1:size(spect, 2), spect_power, [1 size(spect, 2)], [silence_threshold silence_threshold]);
        linkaxes([ax1, ax2], 'x');
    end
end

% figure out silence
not_silent = spect_power >= silence_threshold;

% smooth
p = zeros(1, (silence_smooth_win - 1) / 2);
not_silent_to_smooth = [p double(not_silent) p];
not_silent = conv(not_silent_to_smooth, ones(1, silence_smooth_win), 'valid');
not_silent = not_silent >= silence_smooth_thr;

%% normalize
spect(:, ~not_silent) = 0;
if spect_log
    l = log(spect(:, not_silent));
    l = bsxfun(@minus, l, min(l));
    spect(:, not_silent) = bsxfun(@rdivide, l, sum(l));
else
    spect(:, not_silent) = bsxfun(@rdivide, spect(:, not_silent), sum(spect(:, not_silent)));
end

if debug
    figure;
    imagesc(1:size(spect, 2), 1:size(spect, 1), spect);
end

end

