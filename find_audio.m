function [starts, ends, range_scores] = find_audio(audio, template, fs, varargin)
%FIND_AUDIO

    %% parameters
    alpha = 1.001;
    constrain_length = 0.07; % fraction +/- template length
    match_forward = true;
    match_backward = true;
    fft_window = 512; % samples
    fft_overlap = 472; % samples
    sigma = 2; % ms
    freq_range = [1e3 9e3];
    threshold_score = [];
    debug = false; % popup figures showing scores, path lengths and thresholds
    max_overlap = 0; % fraction of template length

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

    %% prepare audio
    feat_template = features_spectral(template, fs);
    [feat_audio, t_audio] = features_spectral(audio, fs);
    

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
        
        spect = (abs(s(freq_mask, :)) + abs(s2(freq_mask, :))) ./ 2;
    end

    function idx = sliding_window_indices(length, win_length, win_step)
        win = 0:(win_length - 1);
        idx = bsxfun(@plus, 1:win_step:(length - win_length + 1), win');
    end

    function [match_starts, match_ends, match_scores] = match_via_dtw(feat_template, feat_audio, thresh_score)
        % calcualte scores
        [scores, path] = dtw_ua_c(feat_template, feat_audio, alpha);
        
        % debug
        if debug
            o_scores = scores;
        end

        % threshold
        if isempty(thresh_score)
            scores_median = median(scores);
            scores_min = min(scores);
            thresh_score = (scores_median + scores_min) / 2;
        end
        scores(scores >= thresh_score) = nan;

        % threshold by constrain length
        if ~isempty(constrain_length)
            threshold_length = constrain_length * size(feat_template, 2);
            scores(abs(path(1, :) - path(2, :)) > threshold_length) = nan;
        end

        % find minimum by sliding window
        idx = sliding_window_indices(length(scores), round(size(feat_template, 2) / 2), round(size(feat_template, 2) / 4));
        [mn, row] = min(scores(idx));
        time_idx_with_min = row + idx(1, :) - 1; % convert row number to time index

        % potential starts
        potential_ends = unique(time_idx_with_min(~isnan(mn)));
        
        % debug window
        if debug
            % print threshold
            fprintf('Threshold: %f\n', thresh_score);
            
            % plot
            figure;
            ax1 = subplot(4, 1, 1); imagesc(1:size(feat_audio, 2), 1:size(feat_audio, 1), feat_audio);
            ax2 = subplot(4, 1, 2); plot(1:size(feat_audio, 2), o_scores); line([1 size(feat_audio, 2)], [1 1] * thresh_score, 'Color', [1 0 0]);
            ax3 = subplot(4, 1, 3); plot(1:size(feat_audio, 2), abs(path(1, :) - path(2, :))); line([1 size(feat_audio, 2)], [1 1] * threshold_length, 'Color', [1 0 0]);
            ax4 = subplot(4, 1, 4); scatter(find(~isnan(scores)), scores(~isnan(scores)));
            hold on; scatter(potential_ends, scores(potential_ends), 'x'); hold off;
            
            % nice axes
            linkaxes([ax1 ax2 ax3 ax4], 'x');
            h = zoom; h.Motion = 'horizontal'; h.Enable = 'on';
        end

        % corresponding start time
        cor = (1:size(feat_audio, 2)) - size(feat_template, 2) + double(path(1, :));
        
        % prepare return
        match_starts = cor(potential_ends);
        match_ends = potential_ends;
        match_scores = scores(potential_ends);
    end

end