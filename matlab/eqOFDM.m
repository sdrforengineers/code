genOFDMPacket;
% Add random offset
offset = randi([0 1e2]); y = [zeros(offset,1);y]; 
%% Equalization Example
r = awgn(y(offset+1:end),15,'measured'); Fs = 1e6;
pfOffset = comm.PhaseFrequencyOffset('SampleRate',Fs,...
    'FrequencyOffset',Fs*0.001);
r = pfOffset(r);
%% Channel estimation
preambleOFDMMod = comm.OFDMModulator(...
    'FFTLength' ,           FFTLength,...
    'NumGuardBandCarriers', NumGuardBandCarriers,...
    'CyclicPrefixLength',   0,'NumSymbols', 2,'InsertDCNull', true);
od = comm.OFDMDemodulator(preambleOFDMMod);
od.PilotOutputPort = true;
% OFDM Demodulate LLTF
LLTF = Preamble(161:161+160-1); rLLTF = r(161+32:161+160-1);
[rLLTFFreq,rp] = od(rLLTF); [LLTFFreq,p] = od(LLTF(33:end));% remove CP
% Estimate channel
ls = rLLTFFreq./LLTFFreq; % Least-square estimate
chanEst = mean(ls,2); % Average over both symbols
CSI = real(chanEst.*conj(chanEst));
ls = rp./p; % Least-square estimate
chanEstPilots = mean(ls,2); % Average over both symbols
CSIPilots = real(chanEstPilots.*conj(chanEstPilots));
%% Perform Equalization
data = r(2*length(LLTF)+1:end);
odd = comm.OFDMDemodulator(DataOFDMMod);
[dataFreq,pilots] = odd(data);
% Apply LLTF's estimate to data symbols and data pilots
postLLTFEqData = bsxfun(@times, dataFreq, conj(chanEst(:))./CSI(:));
postLLTFEqPilots = ...
    bsxfun(@times, pilots, conj(chanEstPilots(:))./CSIPilots(:));
% Visualization objects
tt1 = comm.ConstellationDiagram;tt2 = comm.ConstellationDiagram;
tt2.Position = tt2.Position + [500 0 0 0];
% Estimate remaining offsets with pilots
correctedSymbols = zeros(size(postLLTFEqData));
for symbol = 1:size(postLLTFEqData,2)
    % Estimate rotation across pilots
    p = postLLTFEqPilots(:,symbol);
    e = conj(mean(p.*conj(Pilots(:,symbol))));
    % Equalize
    sig = postLLTFEqData(:,symbol).*e;
    correctedSymbols(:,symbol) = sig;
    % Visualize
    tt1(sig);tt2(postLLTFEqData(:,symbol));pause(0.1);
end
