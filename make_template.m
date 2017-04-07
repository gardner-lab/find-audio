function [template_fname, source_fname] = make_template()
    % Zebra finch processing
    %
    % Look at audio files (./*.wav). Present a spectrogram, and let the user choose a subset of the spectrogram as an audio template
    % (e.g. for matching using Nathan Perkins's find-audio package (https://github.com/gardner-lab/find-audio)).
    
    fft_size = 512;
    fft_time_shift_seconds = 0.002;
    
    [filename, pathname, ~] = uigetfile('*.wav', 'Audio file containing a suitable template');
    
    source_fname = strcat(pathname, filesep, filename);
    [mic_data, fs] = audioread(source_fname);
    [B, A] = butter(2, [0.03 0.95]);
    mic_data = filtfilt(B, A, double(mic_data));
    
    audio_times = (0:length(mic_data))/fs;
    
    noverlap = fft_size - (floor(fs * fft_time_shift_seconds));
    
    window = hamming(fft_size);
    
    [speck, freqs, times] = spectrogram(mic_data(:,1), window, noverlap, [], fs);
    speck = abs(speck);
    % Adjust "times" to reflect the time at which the information is actually available--i.e. the end,
    % rather than the middle, of the window:
    times = times - times(1) + fft_size/fs;
    
    % Kill off any signal below mean-1std, making the spectrogram prettier:
    speck = reshape(zscore(speck(:)), size(speck)) + 1;
    speck(speck < 0) = 0;
    speck = log(speck + eps);

    f = figure(893476);
    imagesc(times, freqs/1000, speck);
    xlabel('time (s)');
    ylabel('frequency (kHz)');
    axis xy;
    zoom xon;
    
    uicontrol('Style', 'pushbutton', 'String', 'Choose', ...
        'Position', [10 10 80 30], ...
        'BackgroundColor', [0 0.5 0], ...
        'Callback', 'uiresume(gcbf)');
    
    uiwait(gcf);
    
    % If uiwait() returns due to the window getting closed, just exit.
    if ~ishandle(f) || ~isvalid(f)
        disp('Canceled.')
        return;
    end
    
    position = get(gca, 'XLim');
    position_i = audio_times >= position(1) & audio_times <= position(2);
    template = mic_data(position_i);
    
    [filename, pathname, ~] = uiputfile('*_template.wav', 'Save template as...', '0_template.wav');
    template_fname = strcat(pathname, filesep, filename);
    if isequal(filename,0) || isequal(pathname,0)
        disp('Canceled.')
        close(f);
        return;
    end
    
    audiowrite(strcat(pathname, filesep, filename), template, fs);
    fprintf('Saved template as %s\n', strcat(pathname, filesep, filename));
    if ishandle(f) && isvalid(f)
        close(f);
    end
end
