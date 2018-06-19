genOFDMPacket;
% Add random offset with noise
offset = randi([0 1e2]); y = [zeros(offset,1);y];
SNR = 0; r = awgn(y,SNR,'measured');
%% Estimate fine sample offset
LSTF = Preamble(1:160);
LLTF = Preamble(161:161+160-1);
symLen = 80; % FFTLength + CPLen
known = LLTF(1:symLen,1);
% Filter
coeff = conj(flipud(known));
c_filt = abs(filter(coeff, 1, r)).^2;
% Correlation
m = abs(xcorr(r,known)).^2;
padding = length(r)-symLen+1;% Remove padding added to known sequence
c_cor = m(padding:end);
% Plot comparison
stem(c_filt);
hold on; stem(c_cor,'r*');
[v1,peak1] = max(c_cor);
c_cor(peak1) = 0;
[v2,peak2] = max(c_cor);
v = min([v1 v2]);
plot([peak1 peak2],[v v]*0.9,'k','LineWidth',2);hold off;
text(peak1+100,v*0.9,'64 Sample Gap')
grid on;
xlabel('Samples');ylabel('Magnitude');
% Get numerical offset
if abs(peak2-peak1)==FFTLength
    % Adjust position from start of LLTF
    p = min([peak1 peak2]);
    LLTF_Start_est = p-symLen;
    LSTF_Start_est = LLTF_Start_est - length(LSTF);
    % Display
    disp([offset LSTF_Start_est]);
end
