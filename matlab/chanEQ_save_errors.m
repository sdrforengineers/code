h = [0.2; 1; 0.1; 0.02]; % Channel taps
mu = 0.001; % Stepsize
trainingSamples = 1e4;
x = sign(randn(trainingSamples,1)); % Generate BPSK data
r = filter(h,1,x); % Apply channel
L = length(h)+1; h_hat = zeros(L,1);
delta = 3; % Reference tap
%% Estimate channel
index = 1;
[e,x_hat]=deal(zeros(trainingSamples-L+1,1));
for n = L:trainingSamples
    % Select part of training input
    in = r(n:-1:n-L+1);
    % Apply channel estimate to training data
    x_hat(index) = h_hat'*in;
    % Compute error
    e(index) = x(n-delta)-x_hat(index);
    % Update taps
    h_hat = h_hat + mu*conj(e(index))*in;
    index = index + 1;
end
% Slice
x_bar = sign(x_hat);
% Calculate bit errors
eqDelay = L-delta;
fprintf('BER %2.4f\n',mean(x(eqDelay:end-eqDelay-1) ~= x_bar));

%% Plots
figure(1);
semilogy(e.^2);grid on;xlabel('Samples');ylabel('MSE');
figure(2);
Fs = 1e6; N = 64;
htrue=freqz(h,1,N,Fs,'whole');
[hb,we]=freqz(h_hat,1,N,Fs,'whole');
semilogy(we,abs(hb),'b')
hold on;
semilogy(we,abs(htrue),'bo');
semilogy(we,abs(htrue).*abs(hb),'bx');
hold off;
grid on;xlabel('Frequency (Hz)');ylabel('Magnitude');
legend('EQ Response','Channel Response','Composite','Location','Best');


