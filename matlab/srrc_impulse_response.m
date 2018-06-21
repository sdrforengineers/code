close all;
% Use side by side config
config = {'fontsize',16,'fontname','Times'};

% Design filter
betas = [0.1:0.2:0.9];
for beta = betas
span = 16;
sps = 30;
h1 = rcosdesign(beta,span,sps,'sqrt');
N=10;
imp = conv(h1,[zeros(1,N),1,zeros(1,N)]);
il = length(imp);
t = (-il/2:(il-1)/2)./sps;
hold on;plot(t,imp,'DisplayName',['beta=',num2str(beta)]);hold off;
end
legend('show');
%title('Impulse Response of SRRC')
xlabel('t/T_s')
ylabel('Amplitude')
grid on;axis([]);

% Design filter
betas = [0.1:0.2:0.9];
for beta = betas
span = 16;
sps = 30;
h1 = rcosdesign(beta,span,sps,'normal');
N=10;
imp = conv(h1,[zeros(1,N),1,zeros(1,N)]);
il = length(imp);
t = (-il/2:(il-1)/2)./sps;
hold on;plot(t,imp,'DisplayName',['beta=',num2str(beta)]);hold off;
end
legend('show');
%title('Impulse Response of SRRC')
xlabel('t/T_s')
ylabel('Amplitude')
grid on;axis([]);