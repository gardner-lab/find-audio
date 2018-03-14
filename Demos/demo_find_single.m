% Demo of finding single instance of template (uses automatically caulated threshold)
% NOTE: Visualizations require the zftftb library.

%% LOAD

% load audio
[y, fs] = audioread('single.wav');
[template, template_fs] = audioread('template.wav');
if fs ~= template_fs
    error('Mismatched sampling rates.');
end

%% PREPARE VISUAL

% display band
disp_band = [1e3 1e4];

% aligned
[im, f, t] = zftftb_pretty_sonogram(y, fs, 'len', 16.7, 'overlap', 14, ...
    'zeropad', 0, 'filtering', 500, 'clipping', [-2 2], 'norm_amp', 1);

% stard and end index
startidx = find(f<=disp_band(1), 1, 'last');
stopidx = find(f>=disp_band(2), 1, 'first');

% convert to unit8 (0-255)
im = im2uint8(im(startidx:stopidx, :));

% convert to red green blue
im = ind2rgb(im, parula(256));

%% RUN SEARCH

% find
[starts, ends] = find_audio(y, template, fs, 'debug', true, 'match_single', true, 'fft_window', 1024, 'fft_overlap', 1016);
disp([starts ends]);

%% DISPLAY VISUAL

% annotate
figure;
image(t, f, im); axis xy;
line([starts; starts], [ones(size(starts)) * f(1); ones(size(starts)) * f(end)], 'Color', [0 1 0], 'LineWidth', 2);
line([ends; ends], [ones(size(ends)) * f(1); ones(size(ends)) * f(end)], 'Color', [1 0 0], 'LineWidth', 2);

% start zoom
h = zoom;
h.Motion = 'horizontal';
h.Enable = 'on';
