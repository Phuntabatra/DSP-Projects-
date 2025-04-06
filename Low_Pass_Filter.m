% Low Pass Filter for Audio Signals
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

% Normalize audio to [-1, 1] range
audio = audio / max(abs(audio));

%% Filter Parameters
cutoff_freq = 2000; % Cutoff frequency in Hz (adjust as needed)
filter_order = 6;    % Filter order (higher = sharper cutoff)

%% Design Butterworth Low Pass Filter
[b, a] = butter(filter_order, cutoff_freq/(Fs/2), 'low');

%% Filter the Signal
filtered_audio = filtfilt(b, a, audio); % Zero-phase filtering

%% Frequency Response Analysis
freqz(b, a, 1024, Fs);
title(['Frequency Response of ' num2str(cutoff_freq) ' Hz Low Pass Filter']);

%% Time Domain Comparison
t = (0:length(audio)-1)/Fs;
figure;
subplot(2,1,1);
plot(t, audio);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 min(1, length(audio)/Fs)]); % Show first second
grid on;

subplot(2,1,2);
plot(t, filtered_audio);
title(['Filtered Signal (' num2str(cutoff_freq) ' Hz cutoff)']);
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 min(1, length(audio)/Fs)]); % Show first second
grid on;

%% Frequency Domain Comparison
N = length(audio);
f = (0:N/2-1)*Fs/N;

% Original spectrum
audio_fft = abs(fft(audio));
audio_fft = audio_fft(1:N/2);

% Filtered spectrum
filtered_fft = abs(fft(filtered_audio));
filtered_fft = filtered_fft(1:N/2);

figure;
semilogx(f, 20*log10(audio_fft/max(audio_fft)), hold on;
semilogx(f, 20*log10(filtered_fft/max(filtered_fft)), hold off;
title('Frequency Spectrum Comparison');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Original', 'Filtered');
grid on;
xlim([20 Fs/2]);

%% Play Audio
disp('Playing original audio...');
sound(audio, Fs);
pause(length(audio)/Fs + 1); % Wait for playback to finish

disp('Playing filtered audio...');
sound(filtered_audio, Fs);

%% Save Filtered Audio
[~, name, ext] = fileparts(filename);
output_filename = [name '_filtered' ext];
audiowrite(output_filename, filtered_audio, Fs);
disp(['Filtered audio saved as: ' output_filename]);