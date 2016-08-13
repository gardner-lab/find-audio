function threshold = threshold_for_find_audio(varargin)
%THRESHOLD_FOR_FIND_AUDIO Calculate threshold for find_audio
%   This function helps calibrate a threshold for both FIND_AUDIO and
%   FIND_AUDIO_PITCH. This function is meant to be run once for a template
%   and a sample audio clip that includes predominently silence/background
%   noises and a few renditions of the template audio. The function is
%   interactive and does not require any parameters or arguments.
%
%   Note: This function assumes the default parameters for FIND_AUDIO and
%   for FIND_AUDIO_PITCH. If you use custom parameters, you must manually
%   select a threshold. 

%% parameters
template = [];
template_fs = [];
audio = [];
audio_fs = [];
pitch = [];
visualize = true;

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

%% run
% instructions
if isempty(template) || isempty(audio) || isempty(pitch)
    fprintf('The find_audio(...) and find_audio_pitch(...) functions perform most \n');
    fprintf('consistently when they are given a threshold parameter that distinguishes\n');
    fprintf('matches from non-matches. Lower scores are better, so a threshold that is\n');
    fprintf('too low will result in false negatives, and a threshold that is too high\n');
    fprintf('will result in false positives.\n\n');
    fprintf('This function will help calibrate a threshold. You should extract your \n');
    fprintf('template as a WAV file and an example audio file. Ideally, your audio\n');
    fprintf('example should be at least a minute long with representative background\n');
    fprintf('noise and one or more rendition of the template.\n\n');
end

% select template
if isempty(template)
    fprintf('Select template audio file...\n');
    [f, p] = uigetfile({'*.wav;*.ogg;*.flag;*.au;*.aiff;*.aif;*.aifc;*.mp3;*.m4a;*.mp4', 'Audio files'}, 'Pick template audio');
    if isequal(f, 0)
        return;
    end
    [template, template_fs] = audioread(fullfile(p, f));
    fprintf('\n');
end
if isempty(template_fs)
    error('No template sample rate specified.');
end
if 1 < size(template, 2)
    warning('The selected template has more than one channel. This system only supports mono audio and will use the first channel.');
    template(:, 2:end) = [];
end

% select audio
if isempty(audio)
    fprintf('Select sample audio file...\n');
    [f, p] = uigetfile({'*.wav;*.ogg;*.flag;*.au;*.aiff;*.aif;*.aifc;*.mp3;*.m4a;*.mp4', 'Audio files'}, 'Pick sample audio');
    if isequal(f, 0)
        return;
    end
    [audio, audio_fs] = audioread(fullfile(p, f));
    fprintf('\n');
end
if isempty(audio_fs)
    error('No audio sample rate specified.');
end
if 1 < size(audio, 2)
    warning('The selected audio has more than one channel. This system only supports mono audio and will use the first channel.');
    audio(:, 2:end) = [];
end

% check sampling rate
if audio_fs ~= template_fs
    warning('The audio and template have different sample rates. Automatically resampling sample audio.');
    audio = resample(audio, template_fs, audio_fs);
    audio_fs = template_fs;
end

% check relative duration
if length(audio) < 5 * length(template)
    error('The sample audio is too short.');
end
if length(audio) < 50 * length(template)
    warning('Ideally, the sample audio should be much longer than the template. The threshold may not be accurate.');
end

% ask about pitch
if isempty(pitch)
    fprintf('A threshold can be calculated for both the FIND_AUDIO and the FIND_AUDIO_PITCH\n');
    fprintf('functions. FIND_AUDIO is faster, but assumes pretty consistent pitch.\n\n');
    if strcmp(questdlg('Do you expect variability in pitch and plan to use FIND_AUDIO_PITCH? The function is slower, but supports spectral translations.', 'Search Mode', 'No', 'Yes', 'No'), 'Yes')
        pitch = true;
    else
        pitch = false;
    end
end

% calculate threshold
if pitch
    [~, ~, threshold] = find_audio_pitch(audio, template, template_fs, 'calculate_threshold', true);
else
    [~, ~, threshold] = find_audio(audio, template, template_fs, 'calculate_threshold', true);
end

% print threshold
fprintf('Threshold: %.2f\n\n', threshold);

% show matches
if visualize
    fprintf('The figures that are opening will show sample detections based on length and\n');
    fprintf('score of the match, which will give you a sense of the quality of the threshold.\n');
    
    if pitch
        find_audio_pitch(audio, template, template_fs, 'threshold_score', threshold, 'debug', true);
    else
        find_audio(audio, template, template_fs, 'threshold_score', threshold, 'debug', true);
    end
end

end
