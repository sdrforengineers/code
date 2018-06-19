%% General system details
sampleRateHz = 3e6; % Sample rate
samplesPerSymbol = 1;
frameSize = 2^10;
numFrames = 1;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2;
filterUpsample = 4;
filterSymbolSpan = 8;

%% Impairments
snr = 15;
frequencyOffsetHz = 1e5; % Offset in hertz
phaseOffset = 0; % Radians

%% Generate symbols
data = randi([0 samplesPerSymbol], numSamples, 1);
mod = comm.DBPSKModulator();
modulatedData = mod.step(data);

%% Add TX Filter
TxFlt = comm.RaisedCosineTransmitFilter('OutputSamplesPerSymbol', filterUpsample, 'FilterSpanInSymbols', filterSymbolSpan);
filtered = step(TxFlt, modulatedData);

%% Time plot
%figure(1);
%plot(real(modulatedData));
%t = (1:length(filtered))-filterSymbolSpan;
%hold on;plot(t,real(filtered),'r'); hold off;
%xlim([0 300])
%ylim([-2 2])
firInterp = dsp.FIRInterpolator(2);
US = 2;
filtered = resample(filtered,US,1);
%filtered = firInterp(filtered);
modulatedData = resample(modulatedData,US,1);
%modulatedData = firInterp(modulatedData);

df = sampleRateHz/(frameSize*filterUpsample*US);
frequencies2 = -sampleRateHz/2:df:sampleRateHz/2-df;

df = sampleRateHz/(frameSize*US);
frequencies = -sampleRateHz/2:df:sampleRateHz/2-df;

spec = @(sig) fftshift(10*log10(abs(fft(sig))));

%figure(2);
h = plot(frequencies2, spec(filtered),...
     frequencies, spec(modulatedData));
grid on;xlabel('Frequency (Hz)');ylabel('PSD (dB)');
legend('Filtered','Original','Location','Best');
NumTicks = 5;L = h(1).Parent.XLim;
a = sampleRateHz/3;
set(h(1).Parent,'XTick',linspace(-a,a,NumTicks))
xlim([-a a])

