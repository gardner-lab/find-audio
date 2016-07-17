function [warped_time, warped_audio, warped_data] = warp_audio(audio, template, fs, data, varargin)
%WARP_AUDIO


    %% parameters
    alpha = [];
    fft_window = 1024; % samples
    fft_overlap = 1016; % samples
    sigma = 2; % ms
    freq_range = [1e3 9e3];
    debug = false; % popup figures showing scores, path lengths and thresholds

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
    if ~isempty(data) && size(data, 1) ~= length(audio)
        error('Data must be the same sampling rate as the audio.');
    end

    %% prepare audio
    [feat_template, t_template] = features_spectral(template, fs);
    [feat_audio, t_audio] = features_spectral(audio, fs);
    
    %% fill in parameters
    if isempty(alpha)
        cols = size(feat_template, 2);
        alpha = 5 * 0.9 .^ ([0:(ceil(cols / 2) - 1) (floor(cols / 2) - 1):-1:0]) + 2;
    end
    
    %% run
    if debug
        fprintf('Estimated memory usage: %.1f GB\n', (size(feat_audio, 2) + 1) * (size(feat_template, 2) + 1) * 12 / 1024 / 1024 / 1024);
    end
    path = warp_path_via_dtw(feat_template, feat_audio);
    
    % warp timing
    warped_time = [t_template; warp_time(t_audio, path)];
    
    % warp audio
    warped_audio = warp_audio(audio, path);
    warped_audio = exact_rows(warped_audio, length(template), 0);
    
    % warp data
    if isempty(data)
        warped_data = [];
    else
        warped_data = warp_data_by_interp(data, path, length(template));
        %warped_data = warp_data_by_ds(data, path, length(template));
        warped_data = exact_rows(warped_data, length(template), []);
    end
        
    %% helper functions
    
    function [spect, tms] = features_spectral(signal, fs)
        % windows
        t = -fft_window / 2 + 1:fft_window / 2;
        s = fs * sigma / 1e3; % conver sigma from ms to sample count
        w = exp(-(t ./ s).^2); % window
        dw = -2 .* w .* (t ./ (s ^ 2)); % derivative of window
        
        % normalize
        signal = signal ./ max(abs(signal));
        
        % take the two spectrograms, use simple "multi-taper" approach
        [s, freq, tms] = spectrogram(signal, w, fft_overlap, [], fs);
        s2 = spectrogram(signal, dw, fft_overlap, [], fs);
        
        % mask frequency
        freq_mask = freq >= freq_range(1) & freq <= freq_range(2);
        
        spect = (abs(s(freq_mask, :)) + abs(s2(freq_mask, :))) ./ 2;
    end

    function path = warp_path_via_dtw(feat_template, feat_audio)
        % calcualte scores
        path = dtw_path_c(feat_template, feat_audio, alpha);
    end

    function warped = warp_time(time, path)
        warped = zeros(1, path(1, end));
        
        for j = 1:path(1, end)
            idx = path(1, :) == j;
            warped(j) = mean(time(path(2, idx)));
        end
    end

    function warped = warp_audio(signal, path)
        % make indices for FFT
        col = (0:(fft_window - 1))';
        rows = 1:(fft_window - fft_overlap):(1 + length(signal) - fft_window);
        idx = bsxfun(@plus, col, rows);
        
        % making window
        win = hamming(fft_window);
        
        % perform fft
        f = fft(bsxfun(@times, signal(idx), win));
        
        % warped fft
        f_w = zeros(size(f, 1), path(1, end));
        
        for j = 1:path(1, end)
            idx = path(1, :) == j;
            f_w(:, j) = mean(f(:, path(2, idx)), 2);
        end
        
        % inverse fft
        s_w = real(ifft(f_w));
        s_w = bsxfun(@rdivide, s_w, win);
        
        % unwrap
        idx = false(size(s_w));
        idx((fft_overlap / 2):(fft_overlap / 2 + fft_window - fft_overlap - 1), :) = true;
        idx(1:(fft_overlap / 2), 1) = true;
        idx((fft_overlap / 2):end, end) = true;
        
        % produce warped signal
        warped = s_w(idx);
    end

    function warped = warp_data_by_interp(data, path, des_len)
        % get times to sample from (one average signal time per template
        % time)
        smoothed = zeros(1, path(1, end));
        for j = 1:path(1, end)
            idx = path(1, :) == j;
            smoothed(j) = mean(path(2, idx));
        end
        
        % convert smoothed to sample number (from fft sample)
        smoothed = (smoothed - 1) .* (fft_window - fft_overlap) + (fft_window / 2);
        fft_smps = (fft_window / 2) + 8 * (0:(path(1, end) - 1));
        tms = interp1(fft_smps, smoothed, 1:des_len, 'linear', 'extrap');
        
        warped = interp1(1:size(data, 1), data, tms, 'pchip', 'extrap');
        %[warped, ~] = resample(data, tms ./ fs);
    end

    function warped = warp_data_by_ds(data, path, des_len)
        % down sample
        % from # of rows to # of path steps
        data_ds = resample(data, path(2, end), size(data, 1));
        
        % warp data
        data_warped = zeros(path(1, end), size(data, 2));
        for j = 1:path(1, end)
            idx = path(1, :) == j;
            data_warped(j, :) = mean(data_ds(path(2, idx), :), 1);
        end
        
        % up sample
        warped = resample(data_warped, des_len, path(1, end));
    end

    function v = exact_rows(v, rows, pad)
        if size(v, 1) < rows
            if isempty(pad) % use last row
                v = [v; ones(rows - size(v, 1), 1) * v(end, :)];
            else
                v = [v; pad * ones(rows - size(v, 1), size(v, 2), 'like', v)];
            end
        elseif size(v, 1) > rows
            v = v(1:rows, :);
        end
    end
end
