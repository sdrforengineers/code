
qpskMod  = comm.QPSKModulator;
rctFilt = comm.RaisedCosineTransmitFilter( ...
    'OutputSamplesPerSymbol', 2);
varDelay = dsp.VariableFractionalDelay;
awgnChan = comm.AWGNChannel( ...
    'NoiseMethod',  'Signal to noise ratio (SNR)', ...
    'SNR',          25, ...
    'SignalPower',  0.5, ...
    'RandomStream', 'mt19937ar with seed');
rcrFilt = comm.RaisedCosineReceiveFilter( ...
    'InputSamplesPerSymbol',  2, ...
    'DecimationFactor',       1);
rcrFilt2 = comm.RaisedCosineReceiveFilter( ...
    'InputSamplesPerSymbol',  2, ...
    'DecimationFactor',       2);
symsync = comm.SymbolSynchronizer( ...
    'SamplesPerSymbol', 2, ...
    'DampingFactor', sqrt(2)/2, ...
    'NormalizedLoopBandwidth', 0.001);
symsync2 = comm.SymbolSynchronizer( ...
    'SamplesPerSymbol', 2, ...
    'DampingFactor', sqrt(2)/2, ...
    'NormalizedLoopBandwidth', 0.01);


rng(100);
rx = []; rx2 = []; tx = [];
for i = 1:1e3
    data      = randi([0 3], 100, 1);           % Random data
    txSym     = qpskMod(data);                  % QPSK modulation
    txSample  = rctFilt(txSym);                 % Transmit filter
    chanDelay = varDelay(txSample, 0.05*(i-1)); % Variable delay
    chanOut   = awgnChan(chanDelay);            % AWGN channel
    rxSample  = rcrFilt(chanOut);               % Receive filter
    rxSample2  = rcrFilt2(chanOut);               % Receive filter
    rxSym     = symsync(rxSample);              % Symbol synchronizer
    rxSym2     = symsync2(rxSample);              % Symbol synchronizer
    rx = [rx; rxSym];
    rx2 = [rx2; rxSym2];
    tx = [tx;rxSample2];
end

%% Plots
f1 = figure(1);
plot(real(tx),'o');
xlabel('Samples'); ylabel('Amplitude (Real)');grid on;
xlim([0 1e4]);
f2 = figure(2);
plot(real(rx),'.r');
xlabel('Samples'); ylabel('Amplitude (Real)');grid on;
xlim([0 1e4]);
f3 = figure(3);
plot(real(rx2),'xg');
xlim([0 1e4]);
xlabel('Samples'); ylabel('Amplitude (Real)');grid on;

