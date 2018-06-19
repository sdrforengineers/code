% General system details
rng(100);
fs = 1e6; samplesPerSymbol = 1; frameSize = 2^8;
modulationOrder = 2; filterUpsample = 1; filterSymbolSpan = 8;
% Impairments
frequencyOffsetHz = 50000;
% Generate symbols
data = randi([0 samplesPerSymbol], frameSize, 1);
mod = comm.DBPSKModulator(); modulatedData = mod(data);
% Add TX Filter
TxFlt = comm.RaisedCosineTransmitFilter('OutputSamplesPerSymbol',...
    filterUpsample, 'FilterSpanInSymbols', filterSymbolSpan);
filteredData = TxFlt(modulatedData);
filteredData = awgn(filteredData,45,'measured');
% Shift signal in frequency
t = 0:1/fs:(frameSize*filterUpsample-1)/fs;
freqShift = exp(1i.*2*pi*frequencyOffsetHz*t.');   
offsetData = filteredData.*freqShift;
d = offsetData((1:10)+100);
plot(real(d),imag(d),'o');
% Label
for k = 1:length(d)
    text(real(d(k))+0.05,imag(d(k)),num2str(k));
end
hold on;
plot(real(modulatedData),imag(modulatedData),'*r');
hold off; axis([-1 1 -1 1].*1.5)
grid on; ylabel('Quadrature');xlabel('In-phase');
legend('Offset Data','Reference Constellation','Location','NorthEast');

