% General system details
fs = 1e6; samplesPerSymbol = 1; frameSize = 2^8;
modulationOrder = 2; filterUpsample = 4; filterSymbolSpan = 8;
% Impairments
frequencyOffsetHz = 1e5;
% Generate symbols
data = randi([0 samplesPerSymbol], frameSize, 1);
mod = comm.DBPSKModulator(); modulatedData = mod(data);
% Add TX Filter
TxFlt = comm.RaisedCosineTransmitFilter('OutputSamplesPerSymbol',...
    filterUpsample, 'FilterSpanInSymbols', filterSymbolSpan);
filteredData = TxFlt(modulatedData);
filteredData = awgn(filteredData,25,'measured');
% Shift signal in frequency
t = 0:1/fs:(frameSize*filterUpsample-1)/fs;
freqShift = exp(1i.*2*pi*frequencyOffsetHz*t.');   
offsetData = filteredData.*freqShift;
% Estimate
sigNoMod = offsetData.^modulationOrder;
sigNoModOrg = filteredData.^modulationOrder;
% Take FFT and ABS
freqHist = abs(fftshift(fft(sigNoMod)));
freqHistOrg = abs(fftshift(fft(sigNoModOrg)));
% Plot
df = fs/frameSize;
frequencies = -fs/2:df/filterUpsample:fs/2-df/filterUpsample;
spec = @(sig) 10*log10(sig);
% Original
f1 = figure(1);
h = plot(frequencies, spec(freqHistOrg));
grid on;xlabel('Frequency (Hz)');ylabel('PSD^2 (dB)');
NumTicks = 11;L = h(1).Parent.XLim;
set(h(1).Parent,'XTick',linspace(L(1),L(2),NumTicks))
ylim([-10 25]);a = spec(freqHistOrg);hold on;
[~,y] = max(freqHistOrg);
plot(frequencies(y),a(y),'r*','MarkerSize',10);hold off;
% Offset
f2 = figure(2);
h = plot(frequencies, spec(freqHist));
grid on;xlabel('Frequency (Hz)');ylabel('PSD^2 (dB)');
NumTicks = 11;L = h(1).Parent.XLim;
set(h(1).Parent,'XTick',linspace(L(1),L(2),NumTicks))
ylim([-10 25]);a = spec(freqHist);hold on;
[~,y] = max(freqHist);
plot(frequencies(y),a(y),'r*','MarkerSize',10);hold off;
