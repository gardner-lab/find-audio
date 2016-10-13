function [] = make_template();

fft_size = 512;
fft_time_shift_seconds = 0.001;




[filename, pathname, filterindex] = uigetfile('*.wav', 'Audio file containing a suitable template');

[mic_data, fs] = audioread(strcat(pathname, filesep, filename));
[B A] = butter(2, [0.03 0.95]);
mic_data = filtfilt(B, A, double(mic_data));

audio_times = (0:length(mic_data))/fs;

noverlap = fft_size - (floor(fs * fft_time_shift_seconds));

window = hamming(fft_size);

[speck freqs times] = spectrogram(mic_data(:,1), window, noverlap, [], fs);
speck = abs(speck);
% Adjust "times" to reflect the time at which the information is actually available--i.e. the end,
% rather than the middle, of the window:
times = times - times(1) + fft_size/fs;

[nfreqs, ntimes] = size(speck);
med = median(reshape(speck, 1, []));
speck(find(speck <= med)) = med;
f = figure(893475);
im = imagesc(times, freqs/1000, log(speck+eps));
xlabel('time (s)');
ylabel('frequency (kHz)');
axis xy;
zoom xon;

draw_button = uicontrol('Style', 'pushbutton', 'String', 'Choose', ...
    'Position', [10 10 50 20], ...
    'Callback', 'uiresume(gcbf)');
    
uiwait(gcf);

% If uiwait() returns due to the window getting closed, just exit.
if ~ishandle(f) | ~isvalid(f)
    disp('Canceled.')
    return;
end

position = get(gca, 'XLim');
position_i = find(audio_times >= position(1) & audio_times <= position(2));
template = mic_data(position_i);

[filename, pathname, filterindex] = uiputfile('*_template.wav', 'Save template as...', '0_template.wav');
if isequal(filename,0) || isequal(pathname,0)
    disp('Canceled.')
    close(f);
    return;
end

audiowrite(strcat(pathname, filesep, filename), template, fs);
disp(sprintf('Saved template as %s', strcat(pathname, filesep, filename)));
if ishandle(f) & isvalid(f)
    close(f);
end
