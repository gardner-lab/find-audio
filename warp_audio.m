function [warped_time, warped_audio, warped_data, path] = warp_audio(audio, template, fs, data, varargin)
%WARP_AUDIO Uses dynamic time warping to warp audio and data to template
%   This function uses dynamic time warping to find the path that best maps
%   an audio signal onto an audio template, using very similar principles
%   to the FIND_AUDIO command. Spectral features are calculated using the
%   two-taper spectrogram approach with large overlap between windows to
%   provide the best precision during remapping.
%
%   WARPED_TIME = WARP_AUDIO(AUDIO, TEMPLATE, FS) takes an AUDIO clip that
%   was found by FIND_AUDIO (or is in some way already cropped to contain a
%   single song rendition), which is then compared with TEMPLATE through
%   dynamic time warping. Both audio clips must have the same sampling rate
%   FS. WARPED_TIME is a 2 x N matrix where the first row are times
%   throughout the template and the second row represents corresponding
%   times in the new audio signal.
%
%   [WARPED_TIME, WARPED_AUDIO] = WARP_AUDIO(...) generates a new version
%   of AUDIO after warping. Note that this is not playable audio and likely
%   has many artifacts, as audio warping is simply achieved by repeating or
%   removing columns from the STFT.
%
%   [WARPED_TIME, WARPED_AUDIO, WARPED_DATA] = WARP_AUDIO(AUDIO, TEMPLATE, 
%   FS, DATA) will warp auxilary data at the same time was warping the
%   audio. DATA must have the same sampling rate as AUDIO. Warping is
%   achieved by interpolating between points in DATA as needed to add or
%   remove samples. Multiple data streams can be wapred simultaneously, all
%   as columns of DATA.
%
%   [WARPED_TIME, WARPED_AUDIO, WARPED_DATA, PATH] = WARP_AUDIO(...)
%   returns the PATH from the dynamic time wapring process. See DTW_PATH
%   for more details.
%
%   Optional parameters can be passed after the data parameter (if no data
%   to warp, pass an empty matrix for data):
%
%   WARP_AUDIO(..., 'alpha', ALPHA) overrides the default alpha penalty 
%   used for non-diagonal steps in the dynamic time warping. ALPHA can be a
%   scalar or vector corresponding to the number of spectral columns in
%   TEMPLATE. The default values penalizes nondiagonal steps most at the
%   beginning and end of the template.
%
%   WARP_AUDIO(..., 'fft_window', FFT_WINDOW) overides the default size of
%   the window used when calculating spectral features. Defaults to 1024.
%
%   WARP_AUDIO(..., 'fft_overlap', FFT_OVERLAP) overides the default
%   amount of overlap between the window. More overlap provides a more
%   precise warping, but requires more memory and computation. Defaults to
%   1016.
%
%   WARP_AUDIO(..., 'sigma', SIGMA) overides the sigma parameter for the
%   Gaussian window used when calculating the spectral window. Defaults to
%   2ms.
%
%   WARP_AUDIO(..., 'freq_range', FREQ_RANGE) overrides the default
%   frequency range used for warping. FREQ_RANGE must be a vector of length
%   two. Defaults to [1000 9000] Hz.


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
    if nargout >= 2
        warped_audio = warp_audio(audio, path);
        warped_audio = exact_rows(warped_audio, length(template), 0);
    end
    
    % warp data
    if isempty(data) || nargout < 3
        warped_data = [];
    else
        warped_data = warp_data_by_interp(data, path, length(template));
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
        path = dtw_path(feat_template, feat_audio, alpha);
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
        
        % interpolate
        % but leave nan values outside of interpolation range
        % code below replaces that by repeating first/last row
        warped = interp1(1:size(data, 1), data, tms, 'pchip', nan);
        
        % replace nan values
        need_replace = any(isnan(warped), 2);
        
        % add beginning rows
        j = 1;
        row = warped(find(~any(isnan(warped), 2), 1, 'first'), :);
        while need_replace(j)
            warped(j, :) = row;
            j = j + 1;
        end
        
        % add ending rows
        j = size(warped, 1);
        row = warped(find(~any(isnan(warped), 2), 1, 'last'), :);
        while need_replace(j)
            warped(j, :) = row;
            j = j - 1;
        end
        
        %[warped, ~] = resample(data, tms ./ fs);
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
