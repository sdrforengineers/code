clear all;
Fs = 1000;      % Sample rate (Hz)
Fa = 1105;       % Input Frequency (Hz)

% Specify plot parameters
txtsize=10;
ltxtsize=9;
pwidth=4;
pheight=5;
pxoffset=0.65;
pyoffset=0.5;
markersize=5;


% Determine Nyquist zones
zone = 1 + floor(Fa / (Fs/2));
alias = mod(Fa, Fs);
if ~mod(zone,2) % 2nd, 4th, 6th, ... Nyquist Zone
    alias = -(Fs - alias)/Fs;
else            % 3rd, 5th, 7th, ... Nyquist Zone
    alias = (alias)/Fs;
end

N = 2*1/abs(alias) + 1;         % Number of Digital samples
points = 256;                   % Analog points between digital samples
analogIndexes = 0:1/points:N-1;
samplingIndexes = 1:points:length(analogIndexes);
wave = sin(2*pi*Fa/Fs*analogIndexes);
input_str = sprintf('Actual Input Frequency = %i', Fa);
alias_str = sprintf('Measured Frequency = %i', abs(alias * Fs));

% Plot analog signal and sampled
f1=figure(1);
bar(analogIndexes(samplingIndexes), wave(samplingIndexes), .1 ); hold on;
plot(analogIndexes(samplingIndexes), wave(samplingIndexes), 'o')
plot(analogIndexes, wave);
xlabel('Digital Samples');ylabel('Amplitude');
xlim([0 N-1]);
text(1,-.8, input_str,'Color','red','FontSize',12) ;
legend('Digital Sampling', 'Digital Samples', 'Analog Input Signal');
hold off;
ylim([-2 2]);

f2=figure(2)
% Plot analog signal and sampled points
bar(analogIndexes(samplingIndexes), wave(samplingIndexes), .1 ); hold on;
plot(analogIndexes(samplingIndexes), wave(samplingIndexes), 'o')
plot(analogIndexes, sin(2*pi*alias*analogIndexes), '--');

text(1,-.8, alias_str,'Color','red','FontSize',12) ;
legend('Digital Sampling', 'Digital Samples', 'Digital Reconstruction','Location','NorthEast');
xlabel('Digital Samples');ylabel('Amplitude');
xlim([0 N-1])
hold off;
ylim([-2 2]);
