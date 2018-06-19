
rng(105);

%% General system details
sampleRateHz = 1e6; % Sample rate
samplesPerSymbol = 8;
numFrames = 2e2;
modulationOrder = 2;
filterSymbolSpan = 4;
barkerLength = 26; % Must be even

%% Impairments
snr = 20;
% Channel Filter
numTaps = 5;
channelTaps = randn(numTaps,1) + 1i.*randn(numTaps,1);
channelTaps = channelTaps./norm(channelTaps);

%% Generate symbols
bits = randi([0 1],1e4,1);
% Preamble
hBCode = comm.BarkerCode('Length',7,'SamplesPerFrame', barkerLength/2);
barker = step(hBCode)>0;
frame = [barker;barker;bits];
frameSize = length(frame);
modD = comm.QPSKModulator('BitInput',true);
bMod = clone(modD);
modulatedData = 2.*modD.step(frame);

%% Add noise source
chan = comm.AWGNChannel( ...
    'NoiseMethod',  'Signal to noise ratio (SNR)', ...
    'SNR',          snr, ...
    'SignalPower',  3, ...
    'RandomStream', 'mt19937ar with seed');

%% Demodualtor
demod = comm.QPSKDemodulator('BitOutput',true);

%% Equalizer Setup
% Create barker sequence with filtering
barkerMod = bMod.step([barker;barker]);
Mu = 0.00001; % Step size for both estimator and equalizer
numTrainingPoints = length(barkerMod);
% Samples to use for EQ
numTrainingPointsEQ = floor(length(modulatedData)*0.9);
e = zeros(numTrainingPoints-numTaps,1);
indx = 0;indx2 = 0;
e = repmat(e,numFrames,1);
e2 = e;

%% Model of error
BER = zeros(numFrames,1);
PER = zeros(numFrames,1);

preEQ = [];postEQ = [];postEQ2 = [];
prev = 0;
for k=1:numFrames
    
    % Filter with channel filter
    multiPathSignal = filter(channelTaps,1,modulatedData);
    
    % Pass through channel
    noisyData = step(chan, multiPathSignal);
    
    % Pull out preamble
    preambleHat = noisyData(1:length(barkerMod));        
    
   
    %% Equalize Data
    delta = 20; % Tap advantage over channel length (larger is not always better)
    n = numTaps + delta;
    if k==1,f=[1;zeros(n-1,1)];end
    
    for i=n+1:numTrainingPointsEQ
        indx2 = indx2 + 1;
        % Apply equalizer
        y = f'*noisyData(i:-1:i-n+1);
        % Estimate error
        e2(indx2) = modulatedData(i-delta) - y;
        % Update equalizer taps
        f = f + Mu*conj(e2(indx2))*noisyData(i:-1:i-n+1);
    end

    % Apply channel estimate
    frameHatEQ = filter(conj(f),1,noisyData);
    
    % Demodulate and check
    dataHat = demod.step(frameHatEQ(delta+1:end));
    demod.release(); % Reset reference
    errors = abs(dataHat(1:end)-frame(1:length(dataHat)));
    BER(k) = mean(errors);
    PER(k) = BER(k)>0;
    
    % Collect Data
    preEQ = [preEQ; noisyData];
    postEQ = [postEQ; frameHatEQ];

    %% Equalize Data (Blind)
    delta = 20; % Tap advantage over channel length (larger is not always better)
    n = numTaps + delta;
    if k==1,fb=[1;zeros(n-1,1)];end
    
    for i=n+1:numTrainingPointsEQ
        indx2 = indx2 + 1;
        % Apply equalizer
        y = fb'*noisyData(i:-1:i-n+1);
        % Estimate error
        e3 = 2*sign(real(y))+1i*sign(imag(y)) - y;
        % Update equalizer taps
        fb = fb + Mu*conj(e3)*noisyData(i:-1:i-n+1);
    end
    % Apply channel estimate
    frameHatEQ2 = filter(conj(fb),1,noisyData);
    postEQ2 = [postEQ2; frameHatEQ2];
    
end

%% Plot Data Over time
h1 = figure(10);ind = (1:3e2:length(preEQ));
plot(real(preEQ(ind)),'b.');hold on;
plot(imag(preEQ(ind)),'rx');hold off;
grid on;xlabel('Samples');ylabel('Magnitude');legend('I-PreEQ','Q-PreEQ');
ylim([-1 1].*6);
h2 = figure(11);
plot(real(postEQ(ind)),'b.');hold on;
plot(imag(postEQ(ind)),'rx');hold off;
grid on;xlabel('Samples');ylabel('Magnitude');legend('I-PostEQ-LMS','Q-PostEQ-LMS');
ylim([-1 1].*6);
h3 = figure(12);
plot(real(postEQ2(ind)),'b.');hold on;
plot(imag(postEQ2(ind)),'rx');hold off;
grid on;xlabel('Samples');ylabel('Magnitude');legend('I-PostEQ-DD','Q-PostEQ-DD');
ylim([-1 1].*6);

%% Plot loop errors
figure(1);
subplot(3,1,1);
semilogy(abs(e).^2);grid on;xlabel('Samples');ylabel('MSE');title('Chan Error');
subplot(3,1,2);
semilogy(abs(e2).^2);grid on;xlabel('Samples');ylabel('MSE');title('EQ Error');
subplot(3,1,3);
semilogy(BER);grid on;xlabel('Samples');ylabel('BER');title('Frame BER');
% Print some mean results
fprintf('PER %2.5f\n',mean(PER));
fprintf('BER %2.5f\n',mean(BER));

% %% Plot responses
% Fs = 1e6;
% hb=freqz(w,1,1024,sampleRateHz,'whole');
% hf=freqz(f,1,1024,sampleRateHz,'whole');
% [hc,we]=freqz(conv(w,f),1,1024,sampleRateHz,'whole');
% figure(2);
% semilogy(we,abs(hb))
% hold on
% semilogy(we,abs(hf),'r')
% semilogy(we,abs(hc),'go')
% semilogy(we,abs(hb).*abs(hf),'k')
% hold off
% grid on;xlabel('Frequency (Hz)');ylabel('Magnitude');
% %legend('Channel Est','Equalizer Est','Conv(EQ,Chan)','EQ_{PSD} x Chan_{PSD}');
% legend('Channel Est','Equalizer Est','Product Response','Location','Best');





