h = [0.2; 1; 0.1; 0.02]; % Channel taps
mu = 0.001; % Stepsize
trainingSamples = 1e4;
x = sign(randn(trainingSamples,1)); % Generate BPSK data
r = filter(h,1,x); % Apply channel
K = length(h)+2; f = zeros(K,1);
P = length(h)-1; d = zeros(P,1); x_bar_vec = zeros(P,1);
delta = 4; % Reference tap
%% Equalize channel
index = 1;
[e,x_hat]=deal(zeros(trainingSamples-K+1,1));
for n = K:trainingSamples
    % Select part of training input
    in = r(n:-1:n-K+1);
    % Apply channel estimate to training data
    x_hat(index) = f'*in - d'*x_bar_vec;
    % Compute error
    e = x(n-delta)-x_hat(index);
    % Update taps
    f = f + mu*conj(e)*in;
    d = d - mu*conj(e)*x_bar_vec;
    % Update feedback filter
    x_bar_vec = [x(n-delta);x_bar_vec(1:end-1)];
    index = index + 1;
end
% Slice
x_bar = sign(x_hat);
% Calculate bit errors
eqDelay = K-delta;
%fprintf('BER %2.4f\n',mean(x(eqDelay:end-eqDelay-1) ~= x_bar));

%% Check
BERs = [];
for i=1:4
BERs(i) = mean(x(i:i+trainingSamples-K) ~= x_bar);
end
[v,i] = min(BERs)

%% Plots
figure(1);
semilogy(e.^2);grid on;xlabel('Samples');ylabel('MSE');
figure(2);
Fs = 1e6; N = 64;
htrue=freqz(h,1,N,Fs,'whole');
[hb]=freqz(f,1,N,Fs,'whole');
[hd,we]=freqz(d,1,N,Fs,'whole');
semilogy(we,abs(hb),'b')
hold on;
semilogy(we,abs(hd),'b.');
semilogy(we,abs(htrue),'bo');
semilogy(we,abs(htrue).*abs(hb).*abs(hd),'bx');
hold off;
grid on;xlabel('Frequency (Hz)');ylabel('Magnitude');
legend('FF-EQ Response','FB-EQ Response','Channel Response','Composite','Location','Best');


