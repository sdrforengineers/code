% Sinusoid parameters
fs = 1e3; fc = 30; N = 1e3; step = 1;
t = 0:step/fs:(N-1)/fs;
% Create CW Tone
r = cos(2*pi*fc*t); i = sin(2*pi*fc*t);
% Alternatively we can use a hilbert transform from our real signal
y = hilbert(r);
% Estimate frequency from phase
phaseEstHib = unwrap(angle(y))*fs/(2*pi*step); freqEstHib = diff(phaseEstHib);
phaseEstCW = unwrap(atan2(i,r))*fs/(2*pi*step); freqEstCW = diff(phaseEstCW);
tDiff = t(1:end-1);
% Plots
f1 = figure(1);
plot(t,real(y),t,imag(y),'-s',t,r,'*',t,i,'x');
xlabel('Time (Seconds)'); ylabel('Amplitude');grid on;
legend('Real (Hilbert)','Imag (Hilbert)','Cosine','Sine');
xlim([0 60/fs]); ylim([-1.1 1.1])
f2 = figure(2);
plot(tDiff,freqEstHib,tDiff,freqEstCW,'o',...
     tDiff,fc*ones(size(tDiff)),'*');
xlim([0 30/fs]); ylim([0, fc*1.1]);
xlabel('Time (Seconds)');
ylabel('Frequency Estimate (Hz)');
legend('From Hilbert','From Sine/Cosine','True'); grid on;