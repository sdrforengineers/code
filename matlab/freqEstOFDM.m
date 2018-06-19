%% Generate OFDM Waveform
genOFDMPacket;
% Add random offset
offset = randi([0 1e2]); y = [zeros(offset,1);y];
SNR = 20;
%% CFO Estimation
L = 16;  % Short sync field length
m = L;   % Distance between fields
N = 300; % Autocorrelation samples
Fs = 1e6; % Sample rate
% CFO estimation over multiple frequencies
subchannelSpacing = Fs/FFTLength;
CFOs = subchannelSpacing.*(0.01:0.01:0.5);
cfoError = zeros(size(CFOs));
for cfoIndx = 1:length(CFOs)
    % Add noise
    r = awgn(y,SNR,'measured');
    % Frequency offset data
    pfOffset = comm.PhaseFrequencyOffset('SampleRate',Fs,...
        'FrequencyOffset',CFOs(cfoIndx));
    r = pfOffset(r);
    % Determine frequency offsets
    freqEst = zeros(N,1);
    for k=1:N
        P = r(k:k+m)' *  r(k+L:k+m+L);
        freqEst(k) = Fs/L*(angle(P)/(2*pi));
    end
    % Select estimate at offset estimated position
    freqEstimate = freqEst(offset);
    % Calculate CFO error
    cfoError(cfoIndx) = abs(freqEstimate-CFOs(cfoIndx))/subchannelSpacing;
end
% Plot
semilogy(CFOs/subchannelSpacing,cfoError);grid on;
xlabel('CFO Normalized Offset (Hz)');ylabel('CFO Normalized ABS Error (Hz)');
