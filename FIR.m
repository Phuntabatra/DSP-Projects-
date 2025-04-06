% FIR Filter for Audio Signals
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

% If stereo, process each channel separately
if size(audio, 2) == 2
    stereo = true;
    left_ch = audio(:, 1);
    right_ch = audio(:, 2);
else
    stereo = false;
    mono_ch = audio;
end

%% FIR Filter Parameters
filter_type = 'lowpass';  % 'lowpass', 'highpass', 'bandpass', 'bandstop'
cutoff_freq = 2000;       % Cutoff frequency (Hz) for low/highpass
freq_band = [1000 3000];  % Frequency band for bandpass/bandstop (Hz)
filter_order = 101;       % Filter order (number of taps - 1)

%% Design FIR Filter using window method
switch lower(filter_type)
    case 'lowpass'
        fir_coeff = fir1(filter_order, cutoff_freq/(Fs/2), 'low');
    case 'highpass'
        fir_coeff = fir1(filter_order, cutoff_freq/(Fs/2), 'high');
    case 'bandpass'
        fir_coeff = fir1(filter_order, freq_band/(Fs/2), 'bandpass');
    case 'bandstop'
        fir_coeff = fir1(filter_order, freq_band/(Fs/2), 'stop');
    otherwise
        error('Invalid filter type specified');
end

%% Apply FIR Filter
if stereo
    left_filtered = filter(fir_coeff, 1, left_ch);
    right_filtered = filter(fir_coeff, 1, right_ch);
    filtered_audio = [left_filtered, right_filtered];
else
    filtered_audio = filter(fir_coeff, 1, mono_ch);
end

% For zero-phase filtering (uncomment to use):
% if stereo
%     left_filtered = filtfilt(fir_coeff, 1, left_ch);
%     right_filtered = filtfilt(fir_coeff, 1, right_ch);
%     filtered_audio = [left_filtered, right_filtered];
% else
%     filtered_audio = filtfilt(fir_coeff, 1, mono_ch);
% end

%% Time Domain Analysis
t = (0:length(audio)-1)/Fs;
figure('Name', 'FIR Filter Results', 'NumberTitle', 'off');

subplot(3,1,1);
plot(t, audio(:,1), 'b');
title('Original Signal (First Channel)');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 min(0.5, length(audio)/Fs)]); % Show first 0.5 second
grid on;

subplot(3,1,2);
plot(t, filtered_audio(:,1), 'r');
title(['Filtered Signal (' filter_type ', ' num2str(filter_order+1) ' taps)']);
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 min(0.5, length(audio)/Fs)]);
grid on;

% Show impulse response
subplot(3,1,3);
stem(0:filter_order, fir_coeff, 'filled');
title('FIR Filter Coefficients (Impulse Response)');
xlabel('Sample Index');
ylabel('Coefficient Value');
grid on;

%% Frequency Domain Analysis
N = 4096; % Number of points for frequency response
[H, F] = freqz(fir_coeff, 1, N, Fs);

figure;
subplot(2,1,1);
plot(F, 20*log10(abs(H)));
title('FIR Filter Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
xlim([0 Fs/2]);

% Spectrogram comparison
subplot(2,1,2);
window = hamming(512);
noverlap = 256;
nfft = 1024;
spectrogram(filtered_audio(:,1), window, noverlap, nfft, Fs, 'yaxis');
title('Filtered Signal Spectrogram');
colorbar;
ylim([0 8]); % Limit to 8 kHz for typical audio

%% Play Audio
disp('Playing original audio...');
sound(audio, Fs);
pause(length(audio)/Fs + 1);

disp('Playing filtered audio...');
sound(filtered_audio, Fs);

%% Save Filtered Audio
[~, name, ext] = fileparts(filename);
output_filename = [name '_FIR_filtered' ext];
audiowrite(output_filename, filtered_audio, Fs);
disp(['Filtered audio saved as: ' output_filename]);