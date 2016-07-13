% % load audio
% [y, fs] = audioread('audio.wav');
% [template, template_fs] = audioread('template2.wav');
% if fs ~= template_fs
%     error('Mismatched sampling rates.');
% end

% generate audio
fs = 20000;
freq1 = 2093; % C_7
freq2 = 2093 * 2 ^ (2 / 12);
freq3 = 2093 * 2 ^ (4 / 12);
freq4 = 2093 * 2 ^ (6 / 12);
freq5 = 2093 * 2 ^ (8 / 12);
freq6 = 2093 * 2 ^ (10 / 12);
duration = 0.5;
pauses = 1;

t = (0:(1 / fs):duration)';
y = [zeros(pauses * fs, 1); sin(2 * pi * freq1 * t); ...
    zeros(pauses * fs, 1); sin(2 * pi * freq2 * t); ...
    zeros(pauses * fs, 1); sin(2 * pi * freq3 * t); ...
    zeros(pauses * fs, 1); sin(2 * pi * freq4 * t); ...
    zeros(pauses * fs, 1); sin(2 * pi * freq5 * t); ...
    zeros(pauses * fs, 1); sin(2 * pi * freq6 * t); ...
    zeros(pauses * fs, 1)];
y = y + randn(size(y)) / 20;

template = sin(2 * pi * freq1 * t);

% display band
disp_band = [1e3 1e4];

% aligned
[im, f, t] = zftftb_pretty_sonogram(y, fs, ...
    'len', 16.7, 'overlap', 14, 'zeropad', 0, 'filtering', 500, ...
    'clipping', [-2 2], 'norm_amp', 1);

startidx = find(f<=disp_band(1), 1, 'last');
stopidx = find(f>=disp_band(2), 1, 'first');

% convert to unit8 (0-255)
im = im2uint8(im(startidx:stopidx, :));

% convert to red green blue
im = ind2rgb(im, parula(256));

% copy image
im2 = im;

% find
[starts, ends] = find_audio_pitch(y, template, fs, 'debug', true);

% annotate
figure;
image(t, f, im);
line([starts; starts], [ones(size(starts)) * f(1); ones(size(starts)) * f(end)], 'Color', [0 1 0], 'LineWidth', 2);
line([ends; ends], [ones(size(ends)) * f(1); ones(size(ends)) * f(end)], 'Color', [1 0 0], 'LineWidth', 2);

% start zoom
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';