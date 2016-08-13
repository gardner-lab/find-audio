function [starts, ends, range_scores] = find_audio_pitch(audio, template, fs, varargin)
%FIND_AUDIO_PITCH Uses dynamic time warping to find template in audio
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
%   Unlike FIND_AUDIO, this uses reassigned spectrograms with the y-axis
%   log scaled, and allows pitch translations.
%
%   See FIND_AUDIO for usage. This function has the same parameters and
%   behavior.
%
%   There are a few additional parameters avaiable, although they likely do
%   not need to be changed.
%
%   FIND_AUDIO_PITCH(..., 'max_lag', MAX_LAG) sets the maximum vertical
%   translation between the audio and template. Depending on the
%   fundamental frequency, the effect of this varies. The default is 5.
%
%   FIND_AUDIO_PITCH(..., 'log_multiplier', LOG_MULTIPLIER) sets the
%   multiplier on the frequency bin used before taking the log of hte
%   y-axis. The default is 100.

    %% parameters
    alpha = [];
    max_lag = 5;
    log_multiplier = 100;
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
        % normalize
        signal = signal ./ max(abs(signal));
        
        % calculate reassigned spectrogram
        [spect, freq, tms] = ifdv_log(signal, fs, fft_window, fft_overlap, sigma, 1, 1, 5, 5, log_multiplier);
        
        % mask frequency
        freq_mask = freq >= freq_range(1) & freq <= freq_range(2);
        spect = 1 + spect(freq_mask, :);
        
        % log
        if fft_log
            spect = log(spect);
        end
    end

    function idx = sliding_window_indices(length, win_length, win_step)
        win = 0:(win_length - 1);
        idx = bsxfun(@plus, 1:win_step:(length - win_length + 1), win');
    end

    function [match_starts, match_ends, match_scores] = match_via_dtw(feat_template, feat_audio, thresh_score)
        % calcualte scores
        % cor has corresponding start times for minimum warping path
        [scores, cor] = dtpa(feat_template, feat_audio, max_lag, alpha);
        change_in_len = (1:length(scores)) - cor - size(feat_template, 2);
        
        % debug
        if debug
            o_scores = scores;
        end

        % threshold
        if isempty(thresh_score)
            scores_median = median(scores);
            scores_min = min(scores);
            thresh_score = (scores_median + scores_min) / 2;
            
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
        end
        
        % no zero starts
        scores(cor == 0) = nan;

        % find minimum by sliding window
        idx = sliding_window_indices(length(scores), round(size(feat_template, 2) / 2), round(size(feat_template, 2) / 4));
        [mn, row] = min(scores(idx));
        time_idx_with_min = row + idx(1, :) - 1; % convert row number to time index

        % potential starts
        potential_ends = unique(time_idx_with_min(~isnan(mn)));
        
        % debug window
        if debug
            % plot
            figure;
            ax1 = subplot(4, 1, 1); imagesc(1:size(feat_audio, 2), 1:size(feat_audio, 1), feat_audio); axis xy; set(gca, 'Xtick',[], 'Ytick', []);
            ax2 = subplot(4, 1, 2); plot(1:size(feat_audio, 2), o_scores); line([1 size(feat_audio, 2)], [1 1] * thresh_score, 'Color', [1 0 0]); ylabel('Score'); set(gca, 'Xtick',[]);
            ax3 = subplot(4, 1, 3); plot(1:size(feat_audio, 2), abs(change_in_len)); line([1 size(feat_audio, 2)], [1 1] * threshold_length, 'Color', [1 0 0]); ylabel('Score'); set(gca, 'Xtick',[]);
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
