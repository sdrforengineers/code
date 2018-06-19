clear all;
Fs = 1000;      % Sample rate (Hz)
Fa = 1105;      % Input Frequency (Hz)
% Determine Nyquist zones
zone = 1 + floor(Fa / (Fs/2));
alias = mod(Fa, Fs);
if ~mod(zone,2) % 2nd, 4th, 6th, ... Nyquist Zone
    % Its not really a negative amplitude, but it is 180 degrees out of phase,
    % which makes it harder to see on the time domain side, so we cheat
    % to make the graphs look better.
	alias = -(Fs - alias)/Fs;
else            % 3rd, 5th, 7th, ... Nyquist Zone
	alias = (alias)/Fs;
end

% Create the analog/time domain and digital sampling vectors
N = 2*1/abs(alias) + 1;         % Number of Digital samples
points = 256;                   % Analog points between digital samples
analogIndexes = 0:1/points:N-1;
samplingIndexes = 1:points:length(analogIndexes);
wave = sin(2*pi*Fa/Fs*analogIndexes);
% Plot analog  input signal and sampled points
bar(analogIndexes(samplingIndexes), wave(samplingIndexes), .1 ); hold on;
plot(analogIndexes(samplingIndexes), wave(samplingIndexes), 'o')
plot(analogIndexes, wave);
% Plot digitally recreated signal
plot(analogIndexes, sin(2*pi*alias*analogIndexes), '--'); hold off;
legend('Digital Sampling', 'Digital Samples', 'Analog Input', ...
	'Digital Reconstruction', 'Location', 'southoutside');
str = sprintf('Actual Input Frequency = %i\nMeasured Frequency = %i', ...
	Fa, abs(alias * Fs));
title(str);
xlabel('Digital Samples');ylabel('Amplitude');
xlim([0 N-1])
