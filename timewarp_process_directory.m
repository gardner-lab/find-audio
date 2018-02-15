function timewarp_process_directory(template_filename, threshold, minsong, minnonsong)
% Zebra finch song processing:
%
% This is a bridge between raw audio recordings of zebra finch, Nathan Perkins's dynamic timewarp song finder
% (https://github.com/gardner-lab/find-audio), and Ben Pearre's realtime syllable detector
% (https://github.com/gardner-lab/syllable-detector-learn). It reads in the audio, uses dynamic timewarp matching to find songs,
% and saves a song.mat file as required by Ben's learning code.
%
% Apply Nathan Perkins's dynamic timewarp technique to all the wave files in a directory.  I assume they're called
% "channel*.wav".  If the template file '0_template.wav' exists, it will be used, otherwise you will
% be prompted to create it.
%
% template_filename:   filename of a wave file containing one sample of the target song (default 0_template.wav);
%                      if empty, then get the user to define a template.
% threshold:           The match threshold below which stuff is considered a matching song;
%                      the default is to use autothresholding.
% minsong, minnonsong: Keep going until we exceed both of these numbers;
%                      the default is to process all files. If only one is provided, minsong = minnonsong
%
% OUTPUT:              creates a training file, song.mat, as specified in 
%                      https://github.com/gardner-lab/syllable-detector-learn/blob/master/README.md
%
% DEPENDENCIES:        Nathan Perkins's find_audio, https://github.com/gardner-lab/find-audio

AUTO_THRESHOLD_CORRECTION = 1.1; % Bump up the autothreshold by this factor (include more tenuous matches)
NONSONG_THRESHOLD_GAP = 1.2; % Nonsong must be above this factor of threshold
threshold_detect_segment_s = 200; % Build a snippet of audio around this long (seconds) for auto-thresholding
show_detection_points = true;
pause_for_check = false;

files = dir('*.m4a');
[~, sorted_index] = sortrows({files.name}');
% sorted_index = randperm(length(files)); % Maybe randomise the file order?
files = files(sorted_index);
nfiles = length(files);

if ~exist('template_filename', 'var')
    template_filename = '0_template.wav';
end 
if ~exist('minsong', 'var')
    minsong = NaN;
end
if ~exist('minnonsong', 'var')
    if ~isnan(minsong)
        minnonsong = minsong;
    else
        minnonsong = NaN;
    end
end

if isempty(template_filename) ...
        || ( ~exist(template_filename, 'file') && ~exist(strcat(template_filename, '.wav'), 'file'))
    [template_filename, template_source_filename] = make_template();
end
[template, template_fs] = audioread(template_filename);
template_length = length(template);

if ~exist('threshold', 'var') || isempty(threshold) || isnan(threshold) || threshold==0
    randagain = randperm(length(files));

    yn = 'n';
    while strcmp(yn, 'n')
        if exist('template_source_filename', 'var')
            % Include the file from which the template was cut, if that is available
            [template_source, template_source_fs] = audioread(template_source_filename);
        else
            % Otherwise, pick a file randomly; hopefully it will include a match!
            [template_source, template_source_fs] = audioread(files(randagain(1)).name);
        end
        segment_s = length(template_source) / template_source_fs;
        
        i = 1;
        while segment_s < threshold_detect_segment_s
            i = i + 1;
            [a, b] = audioread(files(randagain(i)).name);
            if b ~= template_source_fs
                error('Found an audio file that has different fs than the sample.');
            elseif size(a) == size(template_source) && all(a == template_source)
                continue;
            end
            template_source = [template_source; a];
            segment_s = length(template_source) / template_source_fs;
        end
        figure(412);
        plot_one_spectrogram(template_source, template_source_fs);
        yn = input('Does this spectrogram include a match?  Press y/n <enter>', 's');
    end
    
    
    threshold = threshold_for_find_audio('template', template, 'template_fs', template_fs, 'audio', template_source, 'audio_fs', template_source_fs);
    threshold = threshold * AUTO_THRESHOLD_CORRECTION;
end

song_n = 100;
song_i = 0;
song = zeros(template_length, song_n);

nonsong_n = 100;
nonsong_i = 0;
nonsong = zeros(template_length, nonsong_n);

if ~exist('wbar', 'var') || isempty(wbar) || ~ishandle(wbar) || ~isvalid(wbar)
    wbar = waitbar(0, 'Processing...');
else
    waitbar(0, wbar);
end

start_time = datetime('now');
eta = 'Newtonmas';

for f = 1:nfiles
    if ~isnan(minsong) && ~isnan(minnonsong)
        waitbar(min(song_i/minsong, nonsong_i/minnonsong), wbar, sprintf('Done around %s.', eta));
    else
        waitbar(f/nfiles, wbar, sprintf('Done around %s.', eta));
    end
        
    if song_i >= minsong && nonsong_i >= minnonsong
        % Since the expensive step is find_audio and we have to run that for both song and non-song, there is no point in breaking
        % early if e.g. song_i >> minsong but nonsong_i < minnonsong.
        fprintf('Reached the requested amount of data. %d songs, %d non-songs, %d%% of files processed.\n', ...
            song_i, nonsong_i, round(f*100/nfiles));
        break;
    end
    
    [d, fs] = audioread(files(f).name);
    times = (1:length(d))/fs;
    allsamples_i = [];
    % Find likely matches. The threshold is increased here so we get lots of possible positives, so that the negatives are
    % definitely very much non-song.
    [starts, ends, scores] = find_audio(d, template, fs, 'threshold_score', threshold*NONSONG_THRESHOLD_GAP);
    fprintf('%s: %d regions excluded from non-song, %d above threshold.\n', ...
        files(f).name, length(starts), length(find(scores <= threshold)));
    for n = 1:length(starts)
        % Toss them if they're not good enough.
        % BUT if they're in the grey zone, add them to allsamples_i.
        sample_i = find(times >= starts(n) & times <= ends(n));
        allsamples_i = [allsamples_i sample_i]; % Keep track of these for later...

        if scores(n) <= threshold && song_i < minsong
            sample_length = length(sample_i);
            %% Stretch or compress the audio to make it the same length as the sample? THIS IS A TERRIBLE IDEA!
            %sample_r = resample(sample, template_length, sample_length);
            %% Instead, just pad or chop to fit. This allows the net to learn a few variants.
            sample_under_length = template_length - sample_length;
            if sample_under_length ~= 0 % positive or negative
                begin_fill = ceil(sample_under_length / 2);
                end_fill = floor(sample_under_length / 2);
                sample_i = sample_i(1) - begin_fill : sample_i(end) + end_fill;
            end
            
            % If the padding pushes us over the edge of the available data, discard and ignore.
            if sample_i(1) < 1 || sample_i(end) > length(d)
                continue;
            end
            % Preallocate some more memory, as required
            song_i = song_i + 1;
            if song_i > song_n
                song_n = song_n * 2;
                song(1, song_n) = 0;
            end
            
            song(:, song_i) = d(sample_i);
        end
    end
    
    if nonsong_i <= minnonsong
        non_idx = setdiff(1:length(d), allsamples_i);
        nonsongdata = d(non_idx);
        non_starter = 1;
        while non_starter + template_length - 1 < length(non_idx)
            nonsong_i = nonsong_i + 1;
            if nonsong_i > nonsong_n
                nonsong_n = nonsong_n * 2;
                nonsong(1, nonsong_n) = 0;
            end
            
            % Wrap the nonsong data into the appropriately sized chunks:
            nonsong(:, nonsong_i) = nonsongdata(non_starter:non_starter+template_length-1);
            non_starter = non_starter + template_length;
        end
    end
    
    % Show the detection points
    if show_detection_points
        plot_one_spectrogram(d, fs, threshold, scores, starts, ends);
    end
    
    current_time = datetime('now');
    eta_date = start_time + (current_time - start_time) / min(song_i/minsong, nonsong_i/minnonsong);
    if strcmp(datetime(eta_date, 'Format', 'yyyyMMdd'), datetime(current_time, 'Format', 'yyyyMMdd'))
        eta = datetime(eta_date, 'Format', 'eeee H:mm');
    else
        eta = datetime(eta_date, 'Format', 'H:mm');
    end
    
    if pause_for_check
        fprintf('Press a key...\n');
        pause;
    end
end
close(wbar);

song = song(:, 1:song_i); %#ok<NASGU>
nonsong = nonsong(:, 1:nonsong_i); %#ok<NASGU>

save('song.mat', 'song', 'nonsong', 'fs');



