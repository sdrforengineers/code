clear all;
% Create impulse train of period L and length len with random +/- one
% values
L = 10;
impulse_num = 100; % Total number of impulses in impulse train
len = L*impulse_num;
temp1 = [2*round(rand(impulse_num,1))-1 zeros(impulse_num,L-1)];
x_impulse = reshape(temp1.',[1,L*impulse_num]);
% Create two transmit filter pulse shapes of order L
% Approximate rectangular frequency response --> approximate
% sinc(x) impulse response
txfilt1 = firls(L,[0 0.24 0.25 1],[4 4 0 0]);
% Approximate triangular frequency response --> approximate
% sinc(x)^2 impulse response
txfilt2 = firls(L,[0 0.5 0.52 1],[4 0 0 0]);
% Pulse shape impulse train
y_impulse1 = filter(txfilt1,1,x_impulse);
y_impulse2 = filter(txfilt2,1,x_impulse);

figure(1); plot(x_impulse(L*2:L*4), '-o');
figure(2); plot(y_impulse1(L*2:L*4), '-o');
figure(3); plot(y_impulse2(L*2:L*4), '-o');

% eyediagram(x,n,period,offset)
% creates an eye diagram for the signal x, plotting n samples in each trace
% horizontal axis range between -period/2 and period/2.
eyediagram(y_impulse1,L,L,floor(L/2));
eyediagram(y_impulse2,L,L,floor(L/2));

eyediagram((y_impulse1+0.1*randn(1,length(y_impulse1))),L,L,floor(L/2));
eyediagram((y_impulse2+0.1*randn(1,length(y_impulse2))),L,L,floor(L/2));
