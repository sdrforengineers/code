% Show IQ imbalance

data = randi([0 1], 1e4, 1);
q = comm.QPSKModulator('BitInput',true);
moddata = q(data);

moddata = awgn(moddata,15,'measured');

AmpImb = 2; % 2 dB
PhImb = 25; % 15 degrees
y = iqimbal(moddata,AmpImb,PhImb);
y = y/std(y); % Normalize power


figure(1);
h = scatter(real(moddata),imag(moddata));
grid on;
xlabel('In-phase');ylabel('Quadrature');
h.Parent.FontSize = 14;

figure(2);
h = scatter(real(y),imag(y));
grid on;
xlabel('In-phase');ylabel('Quadrature');
h.Parent.FontSize = 14;
