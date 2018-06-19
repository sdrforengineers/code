% Debugging flags
visuals = true;

%% General system details
sampleRateHz = 1e3; % Sample rate
samplesPerSymbol = 1;
frameSize = 2^10;
numFrames = 30;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2;

%% Setup objects
mod = comm.DBPSKModulator();
cdPre = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Baseband');
cdPost = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'SymbolsToDisplaySource','Property',...
    'SymbolsToDisplay',frameSize/2,...
    'Name','Baseband with Freq Offset');
cdPre.Position(1) = 50;
cdPost.Position(1) = cdPre.Position(1)+cdPre.Position(3)+10;% Place side by side
ap = dsp.ArrayPlot;ap.ShowGrid = true;
ap.Title = 'Frequency Histogram';ap.XLabel = 'Hz';ap.YLabel = 'Magnitude';
ap.XOffset = -sampleRateHz/2;
ap.SampleIncrement = (sampleRateHz)/(2^10);

cdOut = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Baseband');
cdPreOut = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Baseband');

%% Impairments
snr = 15;
frequencyOffsetHz = sampleRateHz*0.02; % Offset in hertz
phaseOffset = 0; % Radians

%% Generate symbols
data = randi([0 samplesPerSymbol], numSamples, 1);
modulatedData = mod.step(data);

%% Add noise
noisyData = awgn(modulatedData,snr);%,'measured');

%% Model of error
% Add frequency offset to baseband signal

% Precalculate constants
normalizedOffset = 1i.*2*pi*frequencyOffsetHz./sampleRateHz;

offsetData = zeros(size(noisyData));
for k=1:frameSize:numSamples
    
    timeIndex = (k:k+frameSize-1).';
    freqShift = exp(normalizedOffset*timeIndex + phaseOffset);
    
    % Offset data and maintain phase between frames
    offsetData(timeIndex) = noisyData(timeIndex).*freqShift;
    
    %     % Visualize Error
    %     if visuals
    %         step(cdPre,noisyData(timeIndex));step(cdPost,offsetData(timeIndex));pause(0.1); %#ok<*UNRCH>
    %     end
    
end

%% Generate Coefficients for PLL
DampingFactor = 1/sqrt(2);
NormalizedLoopBandwidth = 0.01;


% Calculate coefficients based on Michael Rice's analysis
PhaseRecoveryLoopBandwidth = NormalizedLoopBandwidth*samplesPerSymbol;

PhaseRecoveryGain = samplesPerSymbol;

PhaseErrorDetectorGain = 1;

theta = PhaseRecoveryLoopBandwidth/...
    ((DampingFactor + 0.25/DampingFactor)*samplesPerSymbol);

d = 1 + 2*DampingFactor*theta + theta*theta;

% K1
ProportionalGain = (4*DampingFactor*theta/d)/...
    (PhaseErrorDetectorGain*PhaseRecoveryGain);
% K2
IntegratorGain = (4/samplesPerSymbol*theta*theta/d)/...
    (PhaseErrorDetectorGain*PhaseRecoveryGain);

LoopFilter = dsp.IIRFilter( ...
    'Structure', 'Direct form II transposed', ...
    'Numerator', [1 0], 'Denominator', [1 -1]);
Integrator = dsp.IIRFilter(...
    'Structure', 'Direct form II transposed', ...
    'Numerator', [0 1], 'Denominator', [1 -1]);

%% Estimation of error
output = zeros(size(offsetData));
saved_phases = zeros(size(offsetData));
saved_lf = zeros(size(offsetData));

previousSample = complex(0);
loopFiltState = 0;
integFiltState = 0;
DDSPreviousInp = 0;
DigitalSynthesizerGain = -1;
Phase = 0;
FreqEst = 1;

implementation = 'Fast';% 'WithFilters'

for k = 1:length(offsetData)-1
    switch implementation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'Fast'
            % Find Phase Error
            phErr = sign(real(previousSample))*imag(previousSample);
            
            % Phase accumulate and correct
            output(k) = offsetData(k+1)*exp(1i*Phase);
            
            % Loop Filter
            loopFiltOut = phErr*IntegratorGain + loopFiltState;
            loopFiltState = loopFiltOut;
            saved_lf(k) = loopFiltState;
            
            % Direct digital synthesizer implemented as an integrator
            DDSOut = DDSPreviousInp + integFiltState;
            integFiltState = DDSOut;
            DDSPreviousInp = phErr*ProportionalGain+loopFiltOut;
            
            Phase = DigitalSynthesizerGain*DDSOut;
            saved_phases(k) = Phase;
            previousSample = output(k);
            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'WithFilters'
            % Complex phase shift
            output(k) = offsetData(k+1)*exp(1i*Phase);
            
            % Find phase error
            phErr = sign(real(previousSample))*imag(previousSample);
            
            % Loop Filter
            loopFiltOut = step(LoopFilter,phErr*IntegratorGain);
            saved_lf(k) = loopFiltOut;
            
            % Direct Digital Synthesizer
            DDSOut = step(Integrator,phErr*ProportionalGain + loopFiltOut);
            Phase =  DigitalSynthesizerGain * DDSOut;
            saved_phases(k) = Phase;
            previousSample = output(k);
    end
    
end

% %% Static Plot
% H = 1000;k=20000;
% FigHandle = figure;
% set(FigHandle, 'Position', [100, 100, 900, 400]);
% subplot(1,2,1);
% scatter(real(offsetData(k:k+H-1)),imag(offsetData(k:k+H-1)));
% xlim([-1.5 1.5]);ylim([-1.5 1.5]);grid on;xlabel('In-phase');ylabel('Quadrature');
% subplot(1,2,2);
% scatter(real(output(k:k+H-1)),imag(output(k:k+H-1)),'r');
% xlim([-1.5 1.5]);ylim([-1.5 1.5]);grid on;xlabel('In-phase');ylabel('Quadrature');
% 
% %% Animate
% H = 1000;
% for k=1:H:length(output)-H
%     % Visualize Error
%     if visuals
%         step(cdPreOut,offsetData(k:k+H-1));step(cdOut,output(k:k+H-1));pause(0.1); %#ok<*UNRCH>
%     end
% end

%% Determine EVM
symbols = [-1 1];

% Determine closest symbol
est = symbols( (output>0)+1 ).';

% Determine EVM RMS
e = (real(output) - real(est)).^2 + (imag(output) - imag(est)).^2;
top = mean(e); bottom = mean(real(output).^2+imag(output).^2);
evm_rms = sqrt( top/bottom )*100; % Result is in percent
disp(['EVM: ',num2str(evm_rms),'%']);

% %% Static Plots for demonstrating rotating constellation
% r = real(offsetData);
% i = imag(offsetData);
% N = 20;K = 2;
% figure(1);hold on;
% scatter(-1 ,0,'xr');
% scatter(1 ,0,'xr');
% for k=1:K:N
%     scatter(r(k),i(k),'b');
%     text(r(k)+0.02,i(k),num2str(k),'FontSize',12)
% end
% axis([-1.5 1.5 -1.5 1.5])
% grid on;xlabel('In-phase');ylabel('Quadrature');
% legend('Reference Constellation');
% hold off;

% %% Relate phase estimates to loop filter output
% figure;
% plot([diff(saved_phases(1:end-19)),-saved_lf(1:end-20)]);
% grid on;xlabel('Samples');ylabel('Normalized Frequency Offset');
% legend('Diff of Integrator','Loop Filter Output');
% hold off;
%% Determine frequency of offset? (Solution: Used derivative of phase)
% Plot
sig = saved_phases(1:end-100);
x = 1:length(saved_phases);
x2 = 1:length(diff(saved_phases));
%plot(x,saved_phases);
hold on;
plot(-sampleRateHz*diff(sig)/(2*pi),'r');
plot(frequencyOffsetHz.*ones(size(sig)));
hold off;
%ylim([0 frequencyOffsetHz*1.1])
grid on; xlabel('Estimate');ylabel('Offset Hz');






