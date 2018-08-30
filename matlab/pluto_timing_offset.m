% User tunable (samplesPerSymbol>=decimation)
samplesPerSymbol = 12; decimation = 4;
%% System set up
% Set up radio
tx = sdrtx('Pluto','Gain',-10);
rx = sdrrx('Pluto','SamplesPerFrame',1e6,'OutputDataType','double');
% Create binary data for 48, 2-bit symbols
data = randi([0 1],2^15,1);
% Create a QPSK modulator System object with bits as inputs and
qpskMod = comm.QPSKModulator('BitInput',true);
% Modulate and plot the data
modData = qpskMod(data);
% Set up filters
rctFilt = comm.RaisedCosineTransmitFilter( ...
    'OutputSamplesPerSymbol', samplesPerSymbol);
rcrFilt = comm.RaisedCosineReceiveFilter( ...
    'InputSamplesPerSymbol',  samplesPerSymbol, ...
    'DecimationFactor',       decimation);
% Pass data through radio
tx.transmitRepeat(rctFilt(modData));
data = rcrFilt(rx());
% Set up visualization and delay objects
VFD = dsp.VariableFractionalDelay;
cd = comm.ConstellationDiagram;
%% Process received data for timing offset
remainingSPS = samplesPerSymbol/decimation;
% Grab end of data where AGC has converged
data = data(end-remainingSPS*1000+1:end);
for index = 0:300
    % Delay signal
    tau_hat = index/50;
    delayedsig = VFD(data, tau_hat);
    % Linear interpolation
    o = sum(reshape(delayedsig,remainingSPS,...
      length(delayedsig)/remainingSPS).',2)./remainingSPS;
    % Visualize constellation
    cd(o);
    pause(0.1);
end

