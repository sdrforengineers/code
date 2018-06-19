% General system details
rng(100);
fs = 1e6; samplesPerSymbol = 1; frameSize = 2^14;
modulationOrder = 2; filterUpsample = 4; filterSymbolSpan = 8;
% Impairments
frequencyOffsetHz = 50000;
% Generate symbols
data = randi([0 samplesPerSymbol], frameSize, 1);
mod = comm.DBPSKModulator(); modulatedData = mod(data);
% Add TX Filter
TxFlt = comm.RaisedCosineTransmitFilter('OutputSamplesPerSymbol',...
    filterUpsample, 'FilterSpanInSymbols', filterSymbolSpan);
RxFlt = comm.RaisedCosineReceiveFilter('InputSamplesPerSymbol',...
    filterUpsample, 'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', filterUpsample/2);
filteredData = TxFlt(modulatedData);

tx = sdrtx('Pluto','Gain',-30);rx = sdrrx('Pluto',...,
    'OutputDataType','double','SamplesPerFrame',2^16);
rx.CenterFrequency = tx.CenterFrequency + 1e3;

tx.transmitRepeat(filteredData);
%% Capture and sync data
ss = comm.SymbolSynchronizer('TimingErrorDetector','Gardner (non-data-aided)');
for n=1:5
a = RxFlt(rx());
b = ss(a);
end
a = a(end-1024:end);
b = b(end-1024:end);

%% Plot
f1 = figure(1);
plot(real(a),imag(a),'*');
axis([-1 1 -1 1].*1.2)
grid on; ylabel('Quadrature');xlabel('In-phase');
f2 = figure(2);
plot(real(b),imag(b),'*');
axis([-1 1 -1 1].*1.2)
grid on; ylabel('Quadrature');xlabel('In-phase');
%legend('Offset Data','Reference Constellation','Location','NorthEast');
