# DTW audio matching and warping

This repository contains a set of MATLAB tools for matching instances of an auditory template in a longer audio signal, as well as warping the time of two auditory sequences to common timestamps. This is achieved by using varions on dynamic time warping, a dynamic programming technique for mapping two signals together to minimize errors.

In our lab, we use this as a way to extract performances of a zebra finch's song from ongoing recordings. A single template is prepare the screening process, and then the `find_audio` function is able to identify repeated performances of that song quickly and accurately.

This repository offers both MATLAB and C MEX implementations of the actual dynamic programming step, the most computationally expensive part of the process. The MEX files have substantially lower memory usage and are close to a hundred times faster.

## Installation

To install, download the repository and add the main folder to your MATLAB path. You do not need to add subfolders to your path.

The repository currently includes compiled binaries for Mac OS X (x64). If you are using a different platform, you should compile the MEX files for your system to achieve the best performance.

To compile the MEX files, run the following command (assuming you have installed and configured a compiler already):

```
compile_dt_mex
```

## Usage

All usage assumes that you have a **template** audio file, which should be a relatively clean recording of the desired audio recorded under similar circumstances to the signal you will be searching through. (In the case of zebra finches, the template should be a recording of the song with limited cage noises and no socializations recorded from the same cage and microphone where the signal will be coming from.)

### Calibrating the threshold for a template

The `find_audio` and `find_audio_pitch` commands (described below) are both able to automatically determine a threshold, but because that threshold only applies for a single call to the function, it may produce inconsistent results. It is much better to determine a threshold for a template and use that again and again.

The threshold is a cutoff for potential matches. The unanchored dynamic time warping technique used within this repository assigns a score to each potential match. A lower score means the rendition is more similar to the template, and is therefore a better match. Renditions with scores below the threshold are considered, and above the threshold are discarded. Increasing the score may result in more false positives and decreasing the score may result in more false negatives.

It is straight forward to calibrate an initial threshold by providing the template and an audio sample. Ideally, the audio sample should be several times longer than the template and include both several renditions as well as long periods of quiet and background noise.

To determine the threshold, run:

```
threshold_for_find_audio
```

Text written to the command window as well as prompts to select a template and sample audio file will walk you through the process, and print out the identified template.

Note that if you change your template, you should rerun the process above.

### Finding audio

To find instances of a template in an audio stream, use either the command `find_audio` or `find_audio_pitch`. These functions require three parameters, and return three outputs:

```
[starts, ends, scores] = find_audio(audio, template, fs);
```

Parameter `audio` contains the audio stream within which to search for `template`, while `fs` is the sampling rate for both the audio and the template. The function returns `starts` and `ends` times (in seconds) for matches, as well as the `scores` associated with each match.

Assuming you have a threshold for your template as described above, you can pass that into the function as well:

```
[starts, ends, scores] = find_audio(audio, template, fs, 'threshold_score', threshold);
```

There are two different techniques for finding audio. The first (and usually recommended), is `find_audio` which uses dynamic time warping on the two-tapered spectrogram. Given the smoothing in the spectrogram, this will find renditions with small pitch shifts although they will incur a scoring penalty. Alternatively, if you expect more pitch shifting, the `find_audio_pitch` function adds another dimension to the state transition space allowing for translations. It uses the reassigned spectrogram with a log-frequency y-axis, meaning that pitch changes will result in translations of the spectral power.

### Warping audio

In addition to using dynamic time warping to find renditions of the template in an audio signal, this repository includes functionality to use dynamic time warping to warp renditions to match the timing of the template (or just provide equivalent time points between the two).

Unce you have extracted a matching rendition, using the following function to perform the warping:

```
[warped_time, warped_audio, warped_data] = warp_audio(audio, template, fs, data);
```

This funcitons will provide a number of useful outputs regarding the warping needed to match `audio` and `template`.

First, and most broadly useful, `warped_time` is a `2 x N` matrix containing equivalent timestamps from `template` (in the first row) and from `audio` (in the second row). For example, these time stamps can be used to line up experimental data, video frames, or just align audio around a specific syllable or feature.

Second, it provides `warped_audio` which is a warped version of `audio` to match template, although currently this is not playable or useful. The warping is achieved by repeating or removing columns from the spectrogram, which introduces artifacts. That being said, it can be useful for visualizing the audio.

Third, if you provide `data` (with the same sampling rate, `fs`), the function can warp that as well by performing interpolation to add or remove points as needed. The `warped_data` contains this interpolated version.

## Documentation

This readme provides basic information on the usage of the dynamic time warping code, but there are more arguments and parameters available for each of the functions described above. Each function includes extensive documentation in the comments, which can be accessed using the `help` command.

## Known issues

There are a couple of known issues currently:

* The repository only includes MEX binaries for Mac. The MEX files have not yet been tested on other platforms.
* The audio warping is purely cosmetic (it will sound distorted if played). It would be interesting to improve the quality of the warped audio.

## Details

This code is licensed under the MIT license. It was created by [L. Nathan Perkins](https://github.com/nathanntg).

