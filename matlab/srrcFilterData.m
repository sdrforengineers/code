% User tunable (samplesPerSymbol>=decimation)
samplesPerSymbol = 4; decimation = 2;
% Create a QPSK modulator System object and modulate data
qpskMod = comm.QPSKModulator('BitInput',true);
% Set up filters
rctFilt = comm.RaisedCosineTransmitFilter( ...
    'OutputSamplesPerSymbol', samplesPerSymbol);
rcrFilt = comm.RaisedCosineReceiveFilter( ...
    'InputSamplesPerSymbol',  samplesPerSymbol, ...
    'DecimationFactor',       decimation);
% Set up delay object
VFD = dsp.VariableFractionalDelay;
% Delay data with slowly changing delay
rxFilt = [];
for index = 1:1e3
    % Generate, modulate, and tx filter data
    data = randi([0 1],100,1);
    modFiltData = rctFilt(qpskMod(data));
    % Delay signal
    tau_hat = index/30;
    delayedsig = VFD(modFiltData, tau_hat);
    rxSig = awgn(delayedsig,25); % Add noise
    rxFilt = [rxFilt;rcrFilt(rxSig)]; % Rx filter
end