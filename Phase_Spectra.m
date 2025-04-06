% Display Phase Spectra of Audio Signals
clear; clc; close all;

% Load an audio file
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

% Parameters
N = length(audio);
t = (0:N-1)/Fs; % Time vector
f = Fs*(0:(N/2))/N; % Frequency vector

% Compute FFT
audio_fft = fft(audio);
audio_fft = audio_fft(1:N/2+1); % Take single-sided spectrum

% Get magnitude and phase
magnitude = abs(audio_fft);
phase = angle(audio_fft); % Phase in radians

% Unwrap the phase to remove 2*pi jumps
unwrapped_phase = unwrap(phase);

% Plot the original signal
figure;
subplot(3,1,1);
plot(t, audio);
title('Time Domain Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot the magnitude spectrum
subplot(3,1,2);
semilogx(f, 20*log10(magnitude/max(magnitude))); % Plot in dB scale
title('Magnitude Spectrum (dB)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
xlim([20 Fs/2]);

% Plot the phase spectrum
subplot(3,1,3);
plot(f, phase);
title('Phase Spectrum (Wrapped)');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
grid on;
xlim([20 Fs/2]);

% Plot the unwrapped phase spectrum in a separate figure
figure;
plot(f, unwrapped_phase);
title('Unwrapped Phase Spectrum');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
grid on;
xlim([20 Fs/2]);

% Play the original audio
sound(audio, Fs);