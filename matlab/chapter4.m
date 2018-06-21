%----------------------------------------------------------------------
% Chapter 4
% "Digital Communication Systems Engineering Using Software Defined Radio
% MATLAB Scripts
%----------------------------------------------------------------------

% Clear workspace
clear all;

figNum = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sending binary data via sinusoidal signal manipulation

% Parameters
sig_len = 1000; % Signal length (in samples)
sampl_per_bin = 100; % Samples per binary representation
bin_data_len = sig_len/sampl_per_bin; % Length of binary stream is a multiple of signal length
bin_data = round(rand(1,bin_data_len));

% Create sinusoidal carriers
sig_carrier_base = sin(2*pi*(0:(1/sampl_per_bin):(1-(1/sampl_per_bin)))); % Baseline carrier
sig_carrier_freq = sin(2*2*pi*(0:(1/sampl_per_bin):(1-(1/sampl_per_bin)))); % Double frequency
sig_carrier_phase = sin(2*pi*(0:(1/sampl_per_bin):(1-(1/sampl_per_bin)))+(pi/4)); % Phase shifted by 45 degrees

% Modulate sinusoidal carrier via amplitude, phase, and frequency
% manipulations
sig_bin = []; % Binary waveform
sig_ask = []; % Amplitude modulated
sig_psk = []; % Phase modulated
sig_fsk = []; % Frequency modulated
for ind = 1:1:bin_data_len,
    if (bin_data(ind)==1)
        sig_bin = [sig_bin ones(1,sampl_per_bin)];
        sig_ask = [sig_ask sig_carrier_base];
        sig_psk = [sig_psk sig_carrier_base];
        sig_fsk = [sig_fsk sig_carrier_base];
    else
        sig_bin = [sig_bin zeros(1,sampl_per_bin)];
        sig_ask = [sig_ask 0.5*sig_carrier_base];
        sig_psk = [sig_psk sig_carrier_phase];
        sig_fsk = [sig_fsk sig_carrier_freq];
    end;
end;

% Display all three representations
figure(figNum); figNum = figNum+1;
plot(0:1:(sig_len-1),sig_bin);axis([0 sig_len-1 -0.1 1.1]);xlabel('Time (n)');ylabel('Amplitude');
%title('Binary Transmission');

figure(figNum); figNum = figNum+1;
plot(0:1:(sig_len-1),sig_ask);xlabel('Time (n)');ylabel('Amplitude');
%title('Amplitude Shift Keying');

figure(figNum); figNum = figNum+1;
plot(0:1:(sig_len-1),sig_psk);xlabel('Time (n)');ylabel('Amplitude');
%title('Phase Shift Keying');

figure(figNum); figNum = figNum+1;
plot(0:1:(sig_len-1),sig_fsk);xlabel('Time (n)');ylabel('Amplitude');
%title('Frequency Shift Keying');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reducing amount of data transmission using source coding

% Define parameters
len = 10000; % Length of binary data stream

% Have two binary sources, one 50/50 and the other 90/10 in terms of ones
% and zeros
bin1 = round(rand(1,len)); % 50/50 binary
bin2 = round(0.5*rand(1,len)+0.45); %90/10 binary

% Encode strings of ones in terms of the length of these strings
enc_bin1 = [];
enc_bin2 = [];
for ind = 1:1:len,
    if (bin1(ind) == 1)  % Encoding 50/50
        if (ind == 1)
            enc_bin1 = 1;
        else
            enc_bin1(end) = enc_bin1(end)+1;
        end;
    else
        enc_bin1 = [enc_bin1 0];
    end;
    if (bin2(ind) == 1)  % Encoding 90/10
        if (ind == 1)
            enc_bin2 = 1;
        else
            enc_bin2(end) = enc_bin2(end)+1;
        end;
    else
        enc_bin2 = [enc_bin2 0];
    end;
end;

% Find size of encoded binary streams
% (assume all one string length values possess the same number of bits)
ind1 = find(enc_bin1 ~= 0);
ind2 = find(enc_bin2 ~= 0);
[largest_ebin1,ind_largest_ebin1] = max(enc_bin1(ind1));
[largest_ebin2,ind_largest_ebin2] = max(enc_bin2(ind2));
numbits1 = length(dec2bin(largest_ebin1)-'0');
numbits2 = length(dec2bin(largest_ebin2)-'0');
total_size_ebin1 = length(ind1)*numbits1 + length(find(enc_bin1 == 0));
total_size_ebin2 = length(ind2)*numbits2 + length(find(enc_bin2 == 0));

% Plot sizes of original and encoded binary sequences
figure(figNum); figNum = figNum+1;
bar([len total_size_ebin1 total_size_ebin2]);
set(gca,'XTick',[1 2 3],'XTickLabel',{'Original','Encoded (50/50)','Encoded (90/10)'});
xlabel('Binary Transmission Sequences');
ylabel('Number of Bits');
%title('Impact of Source Coding on Binary Transmissions');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Protect data using simple repetition channel coding

% Define parameters
len = 100000; % Length of original binary data stream
N1 = 3; % First repetition factor; should be odd to avoid tie
N2 = 5; % Second repetition factor; should be odd to avoid tie
N3 = 7; % Third repetition factor; should be odd to avoid tie

% Generate binary data stream
bin_str = round(rand(1,len));

% Employ repetition code with repetition factors N1, N2, N3
chcode1_bin_str = zeros(1,N1*len);
chcode2_bin_str = zeros(1,N2*len);
chcode3_bin_str = zeros(1,N3*len);
for ind = 1:1:max([N1 N2 N3]),
    if (ind<=N1)
        chcode1_bin_str(ind:N1:(N1*(len-1)+ind))=bin_str;
    end;
    if (ind<=N2)
        chcode2_bin_str(ind:N2:(N2*(len-1)+ind))=bin_str;
    end;
    if (ind<=N3)
        chcode3_bin_str(ind:N3:(N3*(len-1)+ind))=bin_str;
    end;
end;

% Corrupt both binary strings with zero-mean unit variance Gaussian noise
% followed by rounding (creates "bit flipping" errors)
noisy_bin_str = bin_str + randn(1,len);
rx_bin_str0 = zeros(1,len);
ind0 = find(noisy_bin_str >= 0.5);
rx_bin_str0(ind0) = 1;
noisy_chcode1_bin_str = chcode1_bin_str + randn(1,N1*len);
rx_chcode1_bin_str = zeros(1,N1*len);
ind1 = find(noisy_chcode1_bin_str >= 0.5);
rx_chcode1_bin_str(ind1) = 1;
noisy_chcode2_bin_str = chcode2_bin_str + randn(1,N2*len);
rx_chcode2_bin_str = zeros(1,N2*len);
ind2 = find(noisy_chcode2_bin_str >= 0.5);
rx_chcode2_bin_str(ind2) = 1;
noisy_chcode3_bin_str = chcode3_bin_str + randn(1,N3*len);
rx_chcode3_bin_str = zeros(1,N3*len);
ind3 = find(noisy_chcode3_bin_str >= 0.5);
rx_chcode3_bin_str(ind3) = 1;

% Decode three encoded binary sequences
dec1_bin = (vec2mat(rx_chcode1_bin_str,N1)).';
dec2_bin = (vec2mat(rx_chcode2_bin_str,N2)).';
dec3_bin = (vec2mat(rx_chcode3_bin_str,N3)).';
ind11 = find(((sum(dec1_bin,1))/N1) >= 0.5);
ind12 = find(((sum(dec2_bin,1))/N2) >= 0.5);
ind13 = find(((sum(dec3_bin,1))/N3) >= 0.5);
rx_bin_str1 = zeros(1,len);
rx_bin_str1(ind11) = 1;
rx_bin_str2 = zeros(1,len);
rx_bin_str2(ind12) = 1;
rx_bin_str3 = zeros(1,len);
rx_bin_str3(ind13) = 1;

% Calculate bit error rate
ber0 = sum(abs(bin_str - rx_bin_str0))/len;
ber1 = sum(abs(bin_str - rx_bin_str1))/len;
ber2 = sum(abs(bin_str - rx_bin_str2))/len;
ber3 = sum(abs(bin_str - rx_bin_str3))/len;

% Plot sizes of original and encoded binary sequences
figure(figNum); figNum = figNum+1;

bar([ber0 ber1 ber2 ber3]);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'None',sprintf('N=%d',N1),sprintf('N=%d',N2),sprintf('N=%d',N3)});
xlabel('Binary Transmission Sequences');
ylabel('Bit Error Rate');
%title('Impact of Repetition Coding on Binary Transmissions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating 4-level pulse amplitude modulation, quadrature amplitude modulation,
% and quadrature phase shift keying waveforms

% Define parameters
len = 10000; % Length of binary string
nvar = 0.15; % Noise variance

% Generate the random bit streams that have already been demultiplexed from
% a single high speed data stream
bin_str1 = round(rand(1,len)); % Inphase data stream
bin_str2 = round(rand(1,len)); % Quadrature data stream

% Perform mapping of binary streams to symbols
ind_wavefm = 2.*bin_str2 + 1.*bin_str1; % Create waveform indices
wavefm_4pam = zeros(1,len); % 4-PAM
wavefm_4qam = zeros(1,len); % 4-QAM
wavefm_qpsk = zeros(1,len); % QPSK
symb_4pam = [-3 -1 3 1];
symb_4qam = [-1+i 1+i -1-i 1-i];
symb_qpsk = [exp(i*(pi/5+pi/2)) exp(i*(pi/5+pi)) exp(i*(pi/5+0)) exp(i*(pi/5+3*pi/2)) ];
for ind = 1:1:4,
    wavefm_4pam(find(ind_wavefm == (ind-1))) = symb_4pam(ind);
    wavefm_4qam(find(ind_wavefm == (ind-1))) = symb_4qam(ind);
    wavefm_qpsk(find(ind_wavefm == (ind-1))) = symb_qpsk(ind);
end;

% Add complex zero-mean white Gaussian noise
noise_signal = (1/sqrt(2))*sqrt(nvar)*randn(1,len) + i*(1/sqrt(2))*sqrt(nvar)*randn(1,len);
rx_wavefm_4pam = wavefm_4pam + noise_signal;
rx_wavefm_4qam = wavefm_4qam + noise_signal;
rx_wavefm_qpsk = wavefm_qpsk + noise_signal;

% Generate scatter plots of signal constellations
figure(figNum); figNum = figNum+1;
plot(real(rx_wavefm_4pam),imag(rx_wavefm_4pam),'bo',real(wavefm_4pam),imag(wavefm_4pam),'rx');axis([-4 4 -1 1]);
xlabel('Inphase');ylabel('Quadrature');

figure(figNum); figNum = figNum+1;
plot(real(rx_wavefm_4qam),imag(rx_wavefm_4qam),'bo',real(wavefm_4qam),imag(wavefm_4qam),'rx');axis([-1.5 1.5 -1.5 1.5]);
xlabel('Inphase');ylabel('Quadrature');

figure(figNum); figNum = figNum+1;
plot(real(rx_wavefm_qpsk),imag(rx_wavefm_qpsk),'bo',real(wavefm_qpsk),imag(wavefm_qpsk),'rx');axis([-1.5 1.5 -1.5 1.5]);
xlabel('Inphase');ylabel('Quadrature');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decode messages from previous example using Euclidean distance

% Go through every received waveform and determine Euclidean distance
% between received waveform and the available waveforms
eucl_dist_4pam = zeros(4,len);
eucl_dist_4qam = zeros(4,len);
eucl_dist_qpsk = zeros(4,len);
for ind = 1:1:4,
    eucl_dist_4pam(ind,1:1:len) = abs(symb_4pam(ind).*ones(1,len) - rx_wavefm_4pam);
    eucl_dist_4qam(ind,1:1:len) = abs(symb_4qam(ind).*ones(1,len) - rx_wavefm_4qam);
    eucl_dist_qpsk(ind,1:1:len) = abs(symb_qpsk(ind).*ones(1,len) - rx_wavefm_qpsk);
end;

% Select shortest Euclidean distances
[mdist_4pam,min_ind_4pam] = min(eucl_dist_4pam);
[mdist_4qam,min_ind_4qam] = min(eucl_dist_4qam);
[mdist_qpsk,min_ind_qpsk] = min(eucl_dist_qpsk);

% Decode into estimated binary streams
bin_str_est_4pam = dec2bin(min_ind_4pam-ones(1,len)).';
bin_str_est_4qam = dec2bin(min_ind_4qam-ones(1,len)).';
bin_str_est_qpsk = dec2bin(min_ind_qpsk-ones(1,len)).';

% Calculate bit error rate
ber_4pam = sum([abs((bin_str_est_4pam(1,:)-'0') - bin_str2) ...
	abs((bin_str_est_4pam(2,:)-'0') - bin_str1)])/(2*len);
ber_4qam = sum([abs((bin_str_est_4qam(1,:)-'0') - bin_str2) ...
	abs((bin_str_est_4qam(2,:)-'0') - bin_str1)])/(2*len);
ber_qpsk = sum([abs((bin_str_est_qpsk(1,:)-'0') - bin_str2) ...
	abs((bin_str_est_qpsk(2,:)-'0') - bin_str1)])/(2*len);

% Plot bit error rate results for three modulation schemes
figure(figNum); figNum = figNum+1;
bar([ber_4pam ber_4qam ber_qpsk]);
set(gca,'XTick',[1 2 3],'XTickLabel',{'4-PAM','4-QAM','QPSK'});
xlabel('Modulation Scheme');
ylabel('Bit Error Rate');
%title('Impact of Noise on Modulation Scheme Performance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create "waterfall" curves for the bit error rate of a communication
% system via Monte Carlo

% Define parameters
len = 1000; % Length of individual data transmission
N_snr = 9; % Number of SNR values to evaluation
N_tx = 100; % Number of transmissions per SNR
nvar = [(10.^((1:1:N_snr)/10)).^(-1)]; % Noise variance values

ber_data = zeros(N_snr,N_tx);
for ind = 1:1:N_snr, % Different SNR values
    for ind1 = 1:1:N_tx, % Different transmissions for same SNR value

        % Generate BPSK waveform (we will keep this the same for each
        % SNR value for now)
        tx_sig = 2*round(rand(1,len))-1;

        % Create additive noise
        noise_sig = sqrt(nvar(ind))*randn(1,len);

        % Create received (noisy) signal
        rx_sig = tx_sig + noise_sig;

        % Decode received signal back to binary
        decode_bin_str = zeros(1,len);
        decode_bin_str(find(rx_sig >= 0)) = 1;

        % Determine and store bit error rate
        ber_data(ind,ind1) = sum(abs(decode_bin_str - (tx_sig+1)/2))/len;
    end;
end;

% Calculate mean bit error rate and its standard deviation
mean_ber = mean(ber_data,2).';
std_ber = std(ber_data,'',2).';

% Plot bit error rate waterfall curve ... notice the divergence of standard
% deviation
figure(figNum); figNum = figNum+1;
semilogy(1:1:N_snr,mean_ber,'b-',1:1:N_snr,mean_ber-std_ber,'b--',1:1:N_snr,mean_ber+std_ber,'b-.');
xlabel('Signal to Noise Ratio (dB)');
ylabel('Probability of Error');
%title('Bit Error Rate - BPSK (Mean with Standard Deviation)');
legend('Mean','Mean - Std Dev','Mean + Std Dev');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decoding QAM waveform using I/Q receiver

% Define parameters
N_samp = 1000; % Number of samples per symbol
N_symb = 10; % Number of symbols in transmission
cfreq = 1/10; % Carrier frequency of cosine and sine carriers

% Generate inphase and quadrature channels with 2-PAM waveforms
chI = 2*round(rand(1,N_symb))-1;
chQ = 2*round(rand(1,N_symb))-1;
samp_I = [];
samp_Q = [];
for ind = 1:1:N_symb,
    samp_I = [samp_I chI(ind)*ones(1,N_samp)];
    samp_Q = [samp_Q chQ(ind)*ones(1,N_samp)];
end;

% Apply cosine and sine carriers to inphase and quadrature components, sum
% waveforms together into composite transmission
tx_signal = samp_I.*cos(2.*pi.*cfreq.*(1:1:length(samp_I))) + samp_Q.*sin(2.*pi.*cfreq.*(1:1:length(samp_Q)));

% Separate out inphase and quadrature components from composite
% transmission
sig_I_unfilt = tx_signal.*cos(2.*pi.*cfreq.*(1:1:length(tx_signal)));
sig_Q_unfilt = tx_signal.*sin(2.*pi.*cfreq.*(1:1:length(tx_signal)));
lpf_coeffs = firls(11,[0 cfreq 1.05*cfreq 1],[1 1 0 0]); % Design lowpass filter to remove double frequency term
sig_I_filt = 2.*filter(lpf_coeffs,1,sig_I_unfilt);
sig_Q_filt = 2.*filter(lpf_coeffs,1,sig_Q_unfilt);

% Plot before + after inphase/quadrature signals regarding composite
% waveform creation
figure(figNum); figNum = figNum+1;
plot(0:1:(N_samp*N_symb-1),sig_I_filt,'r',0:1:(N_samp*N_symb-1),samp_I,'b');legend('Recovered','Original');
xlabel('Time (n)');ylabel('Amplitude');
%title('Inphase Signal Component');

figure(figNum); figNum = figNum+1;
plot(0:1:(N_samp*N_symb-1),sig_Q_filt,'r',0:1:(N_samp*N_symb-1),samp_Q,'b');legend('Recovered','Original');
xlabel('Time (n)');ylabel('Amplitude');
%title('Quadrature Signal Component');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gram-Schmidt Orthogonalization and Vectorization

% Define parameters
N_samp = 1000; % Number of samples per time unit

% Create orthonormal basis functions
phi1 = [( 1/sqrt(3))*ones(1,N_samp) ...
	( 1/sqrt(3))*ones(1,N_samp) ...
	(-1/sqrt(3))*ones(1,N_samp)];
phi2 = [( 1/sqrt(6))*ones(1,N_samp) ...
	( 1/sqrt(6))*ones(1,N_samp) ...
	( 2/sqrt(6))*ones(1,N_samp)];
phi3 = [0*ones(1,N_samp) 0*ones(1,N_samp) 0*ones(1,N_samp)];
phi4 = [( 1/sqrt(2))*ones(1,N_samp) ...
	(-1/sqrt(2))*ones(1,N_samp) ...
	  0*ones(1,N_samp)];

% Based on these orthonormal basis functions, create the four symbol waveforms
sig_s1 = (2/sqrt(3))*phi1 + (sqrt(6)/3)*phi2 + 0*phi3 + 0*phi4;
sig_s2 = 0*phi1 + 0*phi2 + 0*phi3 + sqrt(2)*phi4;
sig_s3 = (sqrt(3))*phi1 + 0*phi2 + 0*phi3 + 0*phi4;
sig_s4 = (-1/sqrt(3))*phi1 + (-4/sqrt(6))*phi2 + 0*phi3 + 0*phi4;

% Plot resulting signal waveforms and compare
figure(figNum); figNum = figNum+1;
plot(0:(1/N_samp):(3-(1/N_samp)),sig_s1,'b');xlabel('Time (n)');ylabel('Amplitude');
%title('s_1(n)');
axis([0 3.5 -1.5 1.5]);

figure(figNum); figNum = figNum+1;
plot(0:(1/N_samp):(3-(1/N_samp)),sig_s2,'b');xlabel('Time (n)');ylabel('Amplitude');
%title('s_2(n)');
axis([0 3.5 -1.5 1.5]);

figure(figNum); figNum = figNum+1;
plot(0:(1/N_samp):(3-(1/N_samp)),sig_s3,'b');xlabel('Time (n)');ylabel('Amplitude');
%title('s_3(n)');
axis([0 3.5 -1.5 1.5]);

figure(figNum); figNum = figNum+1;
plot(0:(1/N_samp):(3-(1/N_samp)),sig_s4,'b');xlabel('Time (n)');ylabel('Amplitude');
%title('s_4(n)');
axis([0 3.5 -1.5 1.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlator-based receiver implementation using Gram-Schmidt example
% previous demonstrated

% Define parameters
N_symb = 10; % Number of symbols contained within intercepted signal

% Randomly generate intercepted waveform consisting of s1(n), s2(n), s3(n), and s4(n)
rx_sig = [];
orig_msg = [];
for ind = 1:1:N_symb,
    rnd_val = rand(1,1);
    if (rnd_val < 0.25) % Add s1(n) waveform
        rx_sig = [rx_sig sig_s1];
        orig_msg = [orig_msg 1];
    elseif ((rnd_val >= 0.25)&&(rnd_val < 0.5)) % Add s2(n) waveform
        rx_sig = [rx_sig sig_s2];
        orig_msg = [orig_msg 2];
    elseif ((rnd_val >= 0.5)&&(rnd_val < 0.75)) % Add s3(n) waveform
        rx_sig = [rx_sig sig_s3];
        orig_msg = [orig_msg 3];
    else % Add s4(n) waveform
        rx_sig = [rx_sig sig_s4];
        orig_msg = [orig_msg 4];
    end;
end;

% Vectorize the intercepted signal
dim1_comp = [];
dim2_comp = [];
dim4_comp = [];
for ind = 1:1:N_symb,
    dim1_comp = [dim1_comp sum(rx_sig(((ind-1)*3*N_samp+1):1:(ind*3*N_samp)).*phi1)];
    dim2_comp = [dim2_comp sum(rx_sig(((ind-1)*3*N_samp+1):1:(ind*3*N_samp)).*phi2)];
    dim4_comp = [dim4_comp sum(rx_sig(((ind-1)*3*N_samp+1):1:(ind*3*N_samp)).*phi4)];
end;
dim1_comp = dim1_comp/N_samp;
dim2_comp = dim2_comp/N_samp;
dim4_comp = dim4_comp/N_samp;

% Using the correlator receiver approach, we determine the closest symbol
% vector to each vectorized waveform based on Euclidean distance
s1vec = [(2/sqrt(3)) (sqrt(6)/3) 0 0];
s2vec = [0 0 0 sqrt(2)];
s3vec = [(sqrt(3)) 0 0 0];
s4vec = [(-1/sqrt(3)) (-4/sqrt(6)) 0 0];
est_msg = [];
for ind = 1:1:N_symb,
    [val,symb_ind] = min([ ...
        sum((s1vec - [dim1_comp(ind) dim2_comp(ind) 0 dim4_comp(ind)]).^2) ...
        sum((s2vec - [dim1_comp(ind) dim2_comp(ind) 0 dim4_comp(ind)]).^2) ...
        sum((s3vec - [dim1_comp(ind) dim2_comp(ind) 0 dim4_comp(ind)]).^2) ...
        sum((s4vec - [dim1_comp(ind) dim2_comp(ind) 0 dim4_comp(ind)]).^2) ...
        ]);
    est_msg = [est_msg symb_ind];
end;

% Plot original message and recovered message after intercepted signal was
% vectorized and correlated with available message vectors
figure(figNum); figNum = figNum+1;
stem(0:1:(N_symb-1),orig_msg,'bo');hold on;
stem(0:1:(N_symb-1),est_msg,'rx');
xlabel('Time (n)');ylabel('Message Index');
%title('Matching Correlator Receiver Output with Original Transmission');
legend('Original','Recovered');axis([0 N_symb-1 0 4.5]);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% END OF CHAPTER 4 EXAMPLES %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
