h = [0.5; 1; -0.6]; % Channel to estimate
mu = 0.01; % Stepsize
trainingSamples = 1000;
x = sign(randn(trainingSamples,1)); % Generate BPSK data
r = filter(h,1,x); % Apply channel
L = length(h); h_hat = zeros(L,1);
%% Estimate channel
for n = L:trainingSamples
    % Select part of training input
    in = x(n:-1:n-L+1);
    % Apply channel estimate to training data
    y = h_hat'*in;
    % Compute error
    e = r(n)-y;
    % Update taps
    h_hat = h_hat + mu*conj(e)*in;
end
[h h_hat]

%% Plot responses
Fs = 1e6; N = 64;
htrue=freqz(h,1,N,Fs,'whole');
[hb,we]=freqz(h_hat,1,N,Fs,'whole');
semilogy(we,abs(hb),'b')
hold on;semilogy(we,abs(htrue),'bo');hold off
grid on;xlabel('Frequency (Hz)');ylabel('Magnitude');
legend('Channel Est','Channel Actual','Location','Best');


