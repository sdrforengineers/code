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
% Shift signal in frequency
t = 0:1/fs:(frameSize*filterUpsample-1)/fs;
freqShift = exp(1i.*2*pi*frequencyOffsetHz*t.');   
offsetData = filteredData.*freqShift;
% Plot
df = fs/frameSize;
frequencies = -fs/2:df/filterUpsample:fs/2-df/filterUpsample;
spec = @(sig) fftshift(10*log10(abs(fft(sig))));
% Original
f1 = figure(1);
h = plot(frequencies, spec(filteredData));
grid on;xlabel('Frequency (Hz)');ylabel('PSD (dB)');
NumTicks = 11;L = h(1).Parent.XLim;
set(h(1).Parent,'XTick',linspace(L(1),L(2),NumTicks))
ylim([-20 20]);
% Offset
f2 = figure(2);
h = plot(frequencies, spec(offsetData));
grid on;xlabel('Frequency (Hz)');ylabel('PSD (dB)');
NumTicks = 11;L = h(1).Parent.XLim;
set(h(1).Parent,'XTick',linspace(L(1),L(2),NumTicks))
ylim([-20 20]);
