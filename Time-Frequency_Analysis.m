% Time-Frequency Analysis of Audio Signals
clear; clc; close all;

%% Load Audio File
[filename, pathname] = uigetfile({'*.wav;*.mp3;*.ogg;*.flac;*.m4a;*.mp4',...
                                 'Audio Files (*.wav, *.mp3, *.ogg, *.flac, *.m4a, *.mp4)'},...
                                 'Select an audio file');
if isequal(filename, 0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', fullfile(pathname, filename)]);
end

[audio, Fs] = audioread(fullfile(pathname, filename));

% If stereo, convert to mono by averaging the channels
if size(audio, 2) == 2
    audio = mean(audio, 2);
end

% Normalize audio
audio = audio / max(abs(audio));

%% Parameters
N = length(audio);
t = (0:N-1)/Fs; % Time vector
duration = N/Fs; % Total duration in seconds

%% Time Domain Plot
figure('Name', 'Time-Frequency Analysis', 'NumberTitle', 'off');
subplot(4,1,1);
plot(t, audio);
title('Time Domain Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
xlim([0 duration]);

%% Short-Time Fourier Transform (STFT)
subplot(4,1,2);
window = hamming(1024); % Window function
noverlap = 512; % Overlap between segments
nfft = 1024; % Number of FFT points

spectrogram(audio, window, noverlap, nfft, Fs, 'yaxis');
title('STFT Spectrogram');
colorbar;
ylim([0 8]); % Limit frequency display to 8 kHz for speech/audio

%% Continuous Wavelet Transform (CWT)
subplot(4,1,3);
% For faster computation, we'll analyze a segment if the audio is long
if duration > 5
    segment = audio(1:5*Fs); % First 5 seconds
else
    segment = audio;
end

[cfs, frq] = cwt(segment, Fs);
t_cwt = (0:length(segment)-1)/Fs;
surf(t_cwt, frq, abs(cfs), 'EdgeColor', 'none');
axis tight;
view(0, 90);
title('Continuous Wavelet Transform');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca, 'yscale', 'log');
ylim([20 8000]); % 20 Hz to 8 kHz
colorbar;

