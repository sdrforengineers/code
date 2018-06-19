%% The goal of this file is to demonstrate ISI among OFDM symbols

%% Generate OFDM Waveform
SampleRate = 20e6;
Fs = SampleRate;
PayloadMessage = 'Live long and prosper, from the Communications System Toolbox Team at MathWorks!';
NumFrames = 1;
%% Build OFDM Modulator
FFTLength            = 64;
NumGuardBandCarriers = [6; 5];
NumDataCarriers      = 48;
CyclicPrefixLength   = 16;
PilotCarrierIndices  = [12;26;40;54];
NumOFDMSymInPreamble = 5;
NumBitsPerCharacter  = 7;
% Convert message to bits
msgInBits = double(dec2bin(PayloadMessage, NumBitsPerCharacter).'- 48);
msgInBits = repmat(randi([0 1], NumDataCarriers, 1),10, 1);
PayloadBits = msgInBits(:);
% Calculate number of OFDM symbols per frame
NumOFDMSymbols = ceil(length(PayloadBits)/NumDataCarriers);
% Calculate number of bits padded in each frame
NumPadBits = NumDataCarriers * NumOFDMSymbols - length(PayloadBits);
% Get preamble for each frame
Preamble = double(getOFDMPreambleAndPilot('Preamble', FFTLength, NumGuardBandCarriers));
% Get pilot for each frame
Pilots = double(getOFDMPreambleAndPilot('Pilot', NumOFDMSymbols));
% BPSK modulator
BPSKMod = comm.BPSKModulator;
% OFDM modulator
DataOFDMMod = comm.OFDMModulator(...
    'FFTLength' ,           FFTLength, ...
    'NumGuardBandCarriers', NumGuardBandCarriers, ...
    'InsertDCNull',         true, ...
    'PilotInputPort',       true, ...
    'PilotCarrierIndices',  PilotCarrierIndices, ...
    'CyclicPrefixLength',   CyclicPrefixLength, ...
    'NumSymbols',           NumOFDMSymbols);

%% Modulate Data
symPostBPSK = BPSKMod.step([PayloadBits; randi([0 1], NumPadBits, 1)]);
% OFDM modulation for one frame
symPostOFDM = DataOFDMMod.step(reshape(symPostBPSK, NumDataCarriers, NumOFDMSymbols), Pilots);
% Repeat the frame
y = repmat([Preamble; symPostOFDM], NumFrames, 1);
% Add random offset
offset = randi([0 1e2]);
y = [zeros(offset,1);y];

%% Schmidl and Cox: Coarse Packet Detection
L = 16;  % Short sync field length
m = L;   % Distance between fields
N = 300; % Autocorrelation samples
M = zeros(N,1);
SNR = [0,5,10];
for SNRindx = 1:length(SNR)
    r = awgn(y,SNR(SNRindx),'measured');
    % Determine timing metric
    for k=1:N
        P = r(k:k+m)' *  r(k+L:k+m+L);
        a = abs(y(k+L:k+m+L));
        R = a'*a;
        M(k) = abs(P)^2/(R^2);
    end
    % Plot
    subplot(length(SNR),1,SNRindx);stem(M);
    hold on; stem(offset+1,M(offset+1),'r*'); hold off;
    grid on;xlabel('k');ylabel('M');
    legend('Autocorrelation','True Start');
    title(['SNR: ',num2str(SNR(SNRindx)),'dB']);
end

%% CFO
% L = 16;
% m = L*1;
% N = 300;
% 
% SNR = [0,5,10];
% Fs = 1e6;
% % CFO offset
% subchannelSpacing = Fs/FFTLength;
% CFOs = subchannelSpacing.*(0.01:0.01:0.5);
% cfoError = zeros(size(CFOs));
% figure(2);
% runs = 1e3;
% for cfoIndx = 1:length(CFOs)
%     freqEstimates = zeros(runs,1);
%     pfOffset = comm.PhaseFrequencyOffset('SampleRate',Fs,'FrequencyOffset',CFOs(cfoIndx));
%     parfor run=1:runs
%         r = awgn(y,20,'measured');
%         r = pfOffset(r);
%         % Determine frequency offsets
%         freqEst = zeros(N,1);
%         for k=1:N
%             P = r(k:k+m)' *  r(k+L:k+m+L);
%             freqEst(k) = Fs/L*(angle(P)/(2*pi));
%         end
%         freqEstimates(run) = freqEst(offset);
%     end
%     % Calculate CFO error
%     cfoError(cfoIndx) = abs(mean(freqEstimates)-CFOs(cfoIndx))/subchannelSpacing;
% end
% % Plot
% semilogy(CFOs,cfoError);
% grid on;xlabel('CFO Normalized Offset');ylabel('CFO Normalized ABS Error');
% title(['SNR: ',num2str(SNR(SNRindx)),'dB']);

%% Fine Packet Detection
r = awgn(y,0,'measured');
L = FFTLength;
K = FFTLength/4; % Quarter of short preamble sequence
LLTF = Preamble(161:161+160-1);
LSTF = Preamble(1:160);
symLen = 64;
known = LLTF(1:symLen);

% Filter
coeff = conj(flipud(known));
corr = abs(filter(coeff, 1, r)).^2;

% Correlation
m = abs(xcorr(r,known)).^2;
% Remove padding added to known sequence
padding = length(r)-symLen+1;
m = m(padding:end);

% Plot comparison
figure(3);stem(corr);
hold on; stem(m,'r*');
[v1,peak1] = max(m);
m(peak1) = 0;
[v2,peak2] = max(m);
v = min([v1 v2]);
plot([peak1 peak2],[v v]*0.9,'k','LineWidth',2);hold off;
text(peak1+100,v*0.9,'64 Sample Gap')
grid on;
xlabel('Samples');ylabel('Magnitude');

% Get numerical offset
% [~,peak1] = max(m);
% m(peak1) = 0;
% [~,peak2] = max(m);
if abs(peak2-peak1)==FFTLength
    % Adjust position from start of LLTF
    p = min([peak1 peak2]);
    LLTF_Start_est = p-symLen;
    LSTF_Start_est = LLTF_Start_est - length(LSTF);
    % Display
    disp([offset LSTF_Start_est]);
end

%% Equalization Example
r = awgn(y(LSTF_Start_est+1:end),15,'measured');
pfOffset = comm.PhaseFrequencyOffset('SampleRate',Fs,'FrequencyOffset',1000);
r = pfOffset(r);
%% Channel estimation
preambleOFDMMod = comm.OFDMModulator(...
        'FFTLength' ,           FFTLength,...
        'NumGuardBandCarriers', NumGuardBandCarriers,...
        'CyclicPrefixLength',   0,...
        'NumSymbols', 2,...
        'InsertDCNull', true);
preambleOFDMMod2 = comm.OFDMModulator(...
        'FFTLength' ,           FFTLength,...
        'CyclicPrefixLength',   0,...
        'NumSymbols', 2,...
        'InsertDCNull', true);
od = comm.OFDMDemodulator(preambleOFDMMod);
od.PilotOutputPort = true;
% OFDM Demodulate LLTF
rLLTF = r(161+32:161+160-1);
[rLLTFFreq,rp] = od(rLLTF);
[LLTFFreq,p] = od(LLTF(33:end));% remove CP
% Estimate channel
ls = rLLTFFreq./LLTFFreq; % Least-square estimate   
chanEst = mean(ls,2); % Average over both symbols
CSI = real(chanEst.*conj(chanEst));
ls = rp./p; % Least-square estimate   
chanEstPilots = mean(ls,2); % Average over both symbols
CSIPilots = real(chanEstPilots.*conj(chanEstPilots));
%% Perform Equalization
data = r(321:end);
odd = comm.OFDMDemodulator(DataOFDMMod);
[dataFreq,pilots] = odd(data);
% Apply LLTF's estimate to data symbols and data pilots
postLLTFEqData =  bsxfun(@times, dataFreq, conj(chanEst(:))./CSI(:));
postLLTFEqPilots =  bsxfun(@times, pilots, conj(chanEstPilots(:))./CSIPilots(:));
% Visualization objects
tt1 = comm.ConstellationDiagram;
tt2 = comm.ConstellationDiagram;
% Estimate remaining offsets with pilots
correctedSymbols = zeros(size(postLLTFEqData));
for symbol = 1:size(postLLTFEqData,2)   
    % Estimate rotation across pilots
    p = postLLTFEqPilots(:,symbol);
    x = conj(mean(p.*conj(Pilots(:,symbol))));
    % Equalize
    sig = postLLTFEqData(:,symbol).*x;
    correctedSymbols(:,symbol) = sig;
    % Visualize
    tt1(sig);
    tt2(postLLTFEqData(:,symbol));
    pause(0.1);    
end

%% Demodulate and error check
bd = comm.BPSKDemodulator();
PayloadBitsHat = bd(correctedSymbols(:));
PayloadBitsHatNoEq = bd(postLLTFEqData(:));
PayloadBitsHatNoEqNoPilots = bd(dataFreq(:));

errRate = comm.ErrorRate();
errRateNoEq = comm.ErrorRate();
% Pilots and Eq
errors = errRate(PayloadBits,PayloadBitsHat);
disp(errors(1));
% No Pilots
errors = errRate(PayloadBits,PayloadBitsHatNoEq);
disp(errors(1)); 
% No Eq
errors = errRate(PayloadBits,PayloadBitsHatNoEqNoPilots);
disp(errors(1)); 






