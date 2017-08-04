function [starts, ends, range_scores] = find_audio(audio, template, fs, varargin)
%FIND_AUDIO Uses unanchored dynamic time warping to find template in audio
%   This function uses dynamic time warping to search audio for occurence
%   of a specific template. The approach looks through the audio both
%   forward and backward to find the best possible start and end points,
%   and then applies a two stage thresholding process to find likely
%   candidates. It first uses the variability in the length of the
%   occurence (it must be within a certain percentage of the length of the
%   template). It then looks at the score for the dynamic time warping path
%   as compared with a threshold (the primary tunable parameter for the
%   function).
%
%   After identifying a short list of potential start and end points based
%   on score and length thresholds, the function removes duplicate matches
%   based on overlap, always preserving the match with the best score. 
%
%   [STARTS, ENDS] = FIND_AUDIO(AUDIO, TEMPLATE, FS) searches AUDIO for
%   occurences of TEMPLATE and returns the STARTS and ENDS times (in
%   seconds) of potential matches. Both AUDIO and TEMPLATE must have the
%   same sample rate, FS. STARTS and ENDS are vectors with the same length.
%
%   [STARTS, ENDS, RANGE_SCORES] = FIND_AUDIO(...) also returns the score
%   associated with each match. RANGE_SCORES is a vector of the same length
%   as STARTS and ENDS.
%
%   Optional parameters can be passed as string value pairs. The parameters
%   are (in order of likelihood of use):
%
%   FIND_AUDIO(..., 'threshold_score', THRESHOLD) specifies the threshold
%   for a good match. For consistency, this parameter should always be
%   provided. Helper function THRESHOLD_FOR_FIND_AUDIO helps calculate a
%   good threshold. Dynamic time warping above the threshold are ignored.
%   If not specified, the function will try to calculate a threshold
%   automatically, but this only works well if AUDIO contains both silence
%   and a few renditions of the TEMPLATE audio.
%
%   FIND_AUDIO(..., 'max_overlap', MAX_OVERLAP) determines how much overlap
%   can exist between matches of the template, as a percentage of the
%   template length. The default value is 0.1 (or 10%).
%
%   FIND_AUDIO(..., 'constrain_length', LENGTH) constrains the length of
%   matches to plus or minus the fraction of the template length. LENGTH
%   defaults to 0.07, allowing matches to be plus or minus 7% of the
%   template length.
%
%   FIND_AUDIO(..., 'match_forward', FORWARD) allows you to disable forward
%   matching. By default, the system looks at the audio and template both
%   forward and backwards for the best matches. Using only one direction
%   will be less precise, but faster. Set FORWARD to false to disable.
%
%   FIND_AUDIO(..., 'match_backward', BACKWARD) allows you to disable
%   backward matching. By default, the system looks at the audio and 
%   template both forward and backwards for the best matches. Using only 
%   one direction will be less precise, but faster. Set BACKWARD to false
%   to disable.
%
%   FIND_AUDIO(..., 'debug', DEBUG) set DEBUG to true to enable debugging.
%   Additional information is printed as well as figures opened showing
%   both forward and backward matches.
%
%   FIND_AUDIO(..., 'alpha', ALPHA) overrides the default alpha penalty 
%   used for non-diagonal steps in the dynamic time warping. ALPHA can be a
%   scalar or vector corresponding to the number of spectral columns in
%   TEMPLATE. The default values penalizes nondiagonal steps most at the
%   beginning and end of the template.
%
%   FIND_AUDIO(..., 'fft_window', FFT_WINDOW) overides the default size of
%   the window used when calculating spectral features. Defaults to 512.
%
%   FIND_AUDIO(..., 'fft_overlap', FFT_OVERLAP) overides the default
%   amount of overlap between the window. Defaults to 472.
%
%   FIND_AUDIO(..., 'fft_log', FFT_LOG) allows using the log of the
%   spectrogram for searching. Set FFT_LOG to true to enable.
%
%   FIND_AUDIO(..., 'sigma', SIGMA) overides the sigma parameter for the
%   Gaussian window used when calculating the spectral window. Defaults to
%   2ms.
%
%   FIND_AUDIO(..., 'freq_range', FREQ_RANGE) overrides the default
%   frequency range used for warping. FREQ_RANGE must be a vector of length
%   two. Defaults to [1000 9000] Hz.
%
%   FIND_AUDIO(..., 'match_single', true) conveys that audio contains a
%   single instance of the template. Rather than using a threshold, this
%   causes the matching function to use the lowest start and end match
%   times and will return a single start, end and score.

    %% parameters
    alpha = [];
    constrain_length = 0.07; % fraction +/- template length
    match_forward = true;
    match_backward = true;
    fft_window = 512; % samples
    fft_overlap = 472; % samples
    fft_log = false;
    sigma = 2; % ms
    freq_range = [1e3 9e3];
    threshold_score = [];
    debug = false; % popup figures showing scores, path lengths and thresholds
    max_overlap = 0.1; % fraction of template length
    match_single = false; % audio contains one instance of the template, find it no matter what
    calculate_threshold = false; % used internally for calculating threshold

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
    
    %% prepare return
    starts = [];
    ends = [];
    range_scores = [];
    
    %% checks
    % short circuit a lot of analysis in clear failure cases
    
    % way too short audio
    if length(audio) < length(template) * (1 - constrain_length)
        return;
    end

    %% prepare audio
    feat_template = features_spectral(template, fs);
    [feat_audio, t_audio] = features_spectral(audio, fs);
    
    %% fill in parameters
    if isempty(alpha)
        cols = size(feat_template, 2);
        alpha = 1 * 0.9 .^ ([0:(ceil(cols / 2) - 1) (floor(cols / 2) - 1):-1:0]) + 2;
    end

    %% forward matching
    if match_forward
        [m_ends, m_starts, m_scores] = match_via_dtw(feat_template(:, end:-1:1), feat_audio(:, end:-1:1), threshold_score);
        
        m_ends = 1 + size(feat_audio, 2) - m_ends(end:-1:1);
        m_starts = 1 + size(feat_audio, 2) - m_starts(end:-1:1);
        m_scores = m_scores(end:-1:1);
        
        starts = cat(2, starts, m_starts);
        ends = cat(2, ends, m_ends);
        range_scores = cat(2, range_scores, m_scores);
    end
    
    %% backward matching
    if match_backward
        [m_starts, m_ends, m_scores] = match_via_dtw(feat_template, feat_audio, threshold_score);
        
        starts = cat(2, starts, m_starts);
        ends = cat(2, ends, m_ends);
        range_scores = cat(2, range_scores, m_scores);
    end
    
    %% special mode for calculating thresholds
    if calculate_threshold
        starts = 0;
        ends = 0;
        range_scores = max(range_scores);
        return;
    end
    
    %% deduplicate
    max_overlap_feats = max_overlap * size(feat_template, 2);
    rem = [];
    for k = 1:length(starts)
        overlap = (ends > starts(k) + max_overlap_feats) & (starts < ends(k) - max_overlap_feats);
        if any(range_scores(overlap) < range_scores(k))
            rem = [rem k];
            range_scores(k) = inf;
        end
    end
    starts(rem) = [];
    ends(rem) = [];
    range_scores(rem) = [];
    
    %% sort
    if match_forward && match_backward
        [starts, order] = sort(starts);
        ends = ends(order);
        range_scores = range_scores(order);
    end
    
    %% finalize output
    % convert to times
    starts = t_audio(starts);
    ends = t_audio(ends);

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
        
        % combine spectrogram
        spect = (abs(s(freq_mask, :)) + abs(s2(freq_mask, :))) ./ 2;
        
        % log
        if fft_log
            spect = log(spect);
        end
    end

    function idx = sliding_window_indices(length, win_length, win_step)
        win = 0:(win_length - 1);
        
        % length - win_length + 1 provides all valid indices (but misses
        % some at the end of the length)
        
        % length - win_step + 1 provides all indices (as well as a few past
        % the end of length; as result, clip indices past length)
        
        idx = bsxfun(@plus, 1:win_step:(length - win_step + 1), win');
        
        if idx(end) > length
            idx(idx > length) = length;
        end
    end

    function [match_starts, match_ends, match_scores] = match_via_dtw(feat_template, feat_audio, thresh_score)
        % calcualte scores
        % cor has corresponding start times for minimum warping path
        [scores, cor] = dtw_ua(feat_template, feat_audio, alpha);
        change_in_len = (1:length(scores)) - cor - size(feat_template, 2);
        
        % debug
        if debug
            o_scores = scores;
        end
        
        % single extraction mode
        if match_single
            thresh_score = Inf;
        end

        % threshold
        if isempty(thresh_score)
            scores_median = median(scores);
            scores_min = min(scores);
            thresh_score = 0.75 * scores_median + 0.25 * scores_min;% original= %(scores_median + scores_min) / 2;
            
            if calculate_threshold
                match_starts = 0;
                match_ends = 0;
                match_scores = thresh_score;
                return;
            else
                warning('Threshold determined automatically.');
            end
        end
        scores(scores >= thresh_score) = nan;

        % threshold by constrain length
        if ~isempty(constrain_length)
            threshold_length = constrain_length * size(feat_template, 2);
            scores(abs(change_in_len) > threshold_length) = nan;
        else
            threshold_length = size(feat_template, 2);
        end
        
        if match_single
            % match single... short circuit sliding window logic
            [s, potential_ends] = min(scores);
            
            % no score
            if isnan(s)
                potential_ends = [];
            end
            
            % update start index
            if any(cor(potential_ends) == 0)
                warning('Matches at boundary, may be clipped.');
                cor(cor == 0) = 1;
            end
        else        
            % no zero starts
            scores(cor == 0) = nan;

            % find minimum by sliding window
            idx = sliding_window_indices(length(scores), round(size(feat_template, 2) / 2), round(size(feat_template, 2) / 4));
            [mn, row] = min(scores(idx));
            time_idx_with_min = row + idx(1, :) - 1; % convert row number to time index

            % potential starts
            potential_ends = unique(time_idx_with_min(~isnan(mn)));
        end
        
        % debug window
        if debug
            % plot
            figure;
            ax1 = subplot(4, 1, 1); imagesc(1:size(feat_audio, 2), 1:size(feat_audio, 1), feat_audio); axis xy; set(gca, 'Xtick',[], 'Ytick', []);
            ax2 = subplot(4, 1, 2); plot(1:size(feat_audio, 2), o_scores); line([1 size(feat_audio, 2)], [1 1] * thresh_score, 'Color', [1 0 0]); ylabel('Score'); set(gca, 'Xtick',[]);
            ax3 = subplot(4, 1, 3); plot(1:size(feat_audio, 2), abs(change_in_len)); line([1 size(feat_audio, 2)], [1 1] * threshold_length, 'Color', [1 0 0]); ylabel('Length'); set(gca, 'Xtick',[]);
            ax4 = subplot(4, 1, 4); scatter(find(~isnan(scores)), scores(~isnan(scores))); ylabel('Score'); set(gca, 'Xtick',[]);
            hold on; scatter(potential_ends, scores(potential_ends), 'x'); hold off;
            
            % nice axes
            linkaxes([ax1 ax2 ax3 ax4], 'x');
            h = zoom; h.Motion = 'horizontal'; h.Enable = 'on';
        end
        
        % prepare return
        match_starts = cor(potential_ends);
        match_ends = potential_ends;
        match_scores = scores(potential_ends);
    end

end
