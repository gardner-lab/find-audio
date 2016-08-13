% Demo of warping audio and data (uses automatically caulated threshold)
% and generates simulated data as sin / cos wave.
% NOTE: Visualizations require the zftftb library.

%% LOAD

% load audio
[y, fs] = audioread('template1.wav');
[template, template_fs] = audioread('template2.wav');
if fs ~= template_fs
    error('Mismatched sampling rates.');
end

%% PREPARE VISUAL

% display band
disp_band = [1e3 1e4];

% aligned
[im1, f1, t1] = zftftb_pretty_sonogram(y, fs, ...
    'len', 16.7, 'overlap', 14, 'zeropad', 0, 'filtering', 500, ...
    'clipping', [-2 2], 'norm_amp', 1);

[im2, f2, t2] = zftftb_pretty_sonogram(template, fs, ...
    'len', 16.7, 'overlap', 14, 'zeropad', 0, 'filtering', 500, ...
    'clipping', [-2 2], 'norm_amp', 1);

%% RUN WARPING

% generate some synthetic data
data = [sin((1:length(y)) / 250)' cos((1:length(y)) / 250)'];
[warped_time, warped_audio, warped_data] = warp_audio(y, template, fs, data, 'debug', true);

%% PREPARE VISUAL OF WARPED AUDIO

[im3, f3, t3] = zftftb_pretty_sonogram(warped_audio, fs, ...
    'len', 16.7, 'overlap', 14, 'zeropad', 0, 'filtering', 500, ...
    'clipping', [-2 2], 'norm_amp', 1);

%% DISPLAY VISUAL

figure;
ln = randi(size(warped_time, 2));
ax1 = subplot(3, 1, 1); imagesc(t1, f1, im1); line(warped_time(2, ln) * [1 1], f1([1 length(f1)]), 'Color', [1 1 1]);
ax2 = subplot(3, 1, 2); imagesc(t3, f3, im3); line(warped_time(1, ln) * [1 1], f1([1 length(f1)]), 'Color', [1 1 1]);
ax3 = subplot(3, 1, 3); imagesc(t2, f2, im2);
linkaxes([ax2 ax1 ax3], 'x');

figure;
ax1 = subplot(2, 1, 1); plot((1:length(data)) ./ fs, data); line(warped_time(2, ln) * [1 1], [-1 1], 'Color', [0 1 0]);
ax2 = subplot(2, 1, 2); plot((1:length(warped_data)) ./ fs, warped_data); line(warped_time(1, ln) * [1 1], [-1 1], 'Color', [0 1 0]);
linkaxes([ax2 ax1], 'x');
