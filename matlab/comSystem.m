
%% General system details
sampleRateHz = 1e6; % Sample rate
samplesPerSymbol = 8;
numFrames = 1e3;
modulationOrder = 2;
filterSymbolSpan = 4;
barkerLength = 26; % Must be even

%% Impairments
snr = 15;
% Channel Filter
numTaps = 3;
channelTaps = [1 0.2 -0.1].';
%%channelTaps = [1 0 0].';
%channelTaps = randn(numTaps,1);
channelTaps = channelTaps./max(abs(channelTaps));

%% Generate symbols
bits = randi([0 1],1e3,1); % Generate message (use booktxt.m for a long message)
% Preamble
hBCode = comm.BarkerCode('Length',7,'SamplesPerFrame', barkerLength/2);
barker = step(hBCode)>0;
frame = [barker;barker;bits];
frameSize = length(frame);
modD = comm.DBPSKModulator();
bMod = clone(modD);
modulatedData = modD.step(frame);

%% Add TX/RX Filters
TxFlt = comm.RaisedCosineTransmitFilter(...
    'OutputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan);

RxFlt = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', samplesPerSymbol);% Set to filterUpsample/2 when introducing timing estimation
TxFltRef = clone(TxFlt);
RxFltRef = clone(RxFlt);

%% Add noise source
chan = comm.AWGNChannel( ...
    'NoiseMethod',  'Signal to noise ratio (SNR)', ...
    'SNR',          snr, ...
    'SignalPower',  1, ...
    'RandomStream', 'mt19937ar with seed');

%% Demodualtor
demod = comm.DBPSKDemodulator;

%% Equalizer Setup
% Create barker sequence with filtering
barkerMod = bMod.step([barker;barker]);
sig = [TxFltRef.step(barkerMod);...
    zeros(TxFltRef.OutputSamplesPerSymbol*TxFltRef.FilterSpanInSymbols,1)];
% Pass through filters and remove lag
barkerRef = RxFltRef.step(sig);
barkerRef = barkerRef(filterSymbolSpan+1:...
                      filterSymbolSpan+length(barkerMod));

Mu = 0.01; % Step size for both estimator and equalizer

numTrainingPoints = length(barkerRef);
e = zeros(numTrainingPoints-numTaps,1);
indx = 0;indx2 = 0;
e = repmat(e,numFrames,1);
e2 = e;

%% Setup visualization object(s)
hAP = dsp.ArrayPlot;
hAP.NumInputPorts = 2;
mb = min(channelTaps) - abs(min(channelTaps))*0.5;
mt = max(channelTaps) + abs(max(channelTaps))*0.5;
hAP.YLimits = [mb mt];

%% Model of error
BER = zeros(numFrames,1);
PER = zeros(numFrames,1);

for k=1:numFrames
    
    % Pad to get full signal from filter
    modulatedDataPad = [modulatedData; zeros(filterSymbolSpan*2,1)];
    
    % Filter signal
    filteredTXDataDelayed = step(TxFlt, modulatedDataPad);
    
    % Filter with channel filter
    multiPathSignal = filter(upsample(channelTaps,samplesPerSymbol)...
        ,1,filteredTXDataDelayed);
    
    % Add AWGN
    noisyData = step(chan, multiPathSignal);
    
    % RxFilter signal
    filteredData = step(RxFlt, noisyData);
    
    % Pull out preamble
    delay = filterSymbolSpan + 1;
    preambleHat = filteredData(delay:delay+length(barkerRef));        
        
    %% Equalize Data (Based on Preamble)
    delta = 1; % Tap advantage over channel length (larger is not always better)
    n = numTaps + delta;
    if k==1,f=zeros(1,n);end
    
    for i=n+1:numTrainingPoints
        indx2 = indx2 + 1;
        % Apply equalizer
        y = f*preambleHat(i:-1:i-n+1);
        % Estimate error
        e2(indx2) = barkerMod(i-delta) - y;
        % Update equalizer taps
        f = f + Mu*e2(indx2)*preambleHat(i:-1:i-n+1)';
    end
    
    % Extract frame
    frameHat = filteredData(delay:delay + length(modulatedData) - 1);
    
    % Apply channel estimate
    frameHatEQ = filter(f,1,frameHat);
   
    % Demodulate and check
    dataHat = demod.step(frameHatEQ(delta+1:end));
    demod.release(); % Reset reference
    %[dataHat frame(1:length(dataHat))]
    errors = abs(dataHat-frame(1:length(dataHat)));
    BER(k) = mean(errors);
    PER(k) = BER(k)>0;
end

%% Plot loop errors
figure(1);
%subplot(3,1,1);
%semilogy(abs(e));grid on;xlabel('Samples');ylabel('RMSE');title('Chan Error');
%subplot(3,1,2);
semilogy(abs(e2));grid on;xlabel('Samples');ylabel('RMSE');title('EQ Error');
%subplot(3,1,3);
%semilogy(BER);grid on;xlabel('Samples');ylabel('BER');title('Frame BER');
% Print some mean results
fprintf('PER %2.5f\n',mean(PER));
fprintf('BER %2.5f\n',mean(BER));

%% Plot responses
Fs = 1e6; N = 64;
htrue=freqz(channelTaps,1,N,sampleRateHz,'whole');
%hb=freqz(w,1,N,sampleRateHz,'whole');
[hf,we]=freqz(f,1,N,sampleRateHz,'whole');
%[hc,we]=freqz(conv(w,f),1,N,sampleRateHz,'whole');
figure(2);
%semilogy(we,abs(hb),'b')
hold on
semilogy(we,abs(htrue),'bo')
semilogy(we,abs(hf),'r')
%semilogy(we,abs(hc),'go')
%semilogy(we,abs(hb).*abs(hf),'k')
semilogy(we,abs(htrue).*abs(hf),'k*')
hold off
grid on;xlabel('Frequency (Hz)');ylabel('Magnitude');
legend('Channel Est','Channel Actual','Equalizer Est','Combined','Combined Actual','Location','Best');
%legend('Channel Est','Equalizer Est','Product Response','Location','Best');





