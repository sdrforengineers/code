% Create figure to describe OFDM symbol structure
nFFT = 64;
subcarrierIndexes = 1:nFFT;
datacarrierIndexes = 8:64-7;

pilotTones = [11,25,40,54];
carriers = zeros(nFFT,1);
guardcarriers = subcarrierIndexes;
guardcarriers(datacarrierIndexes) = [];


% Plot
stem(guardcarriers-1, carriers(guardcarriers),'rx');
hold on;
carriers(datacarrierIndexes) = 1;
stem(datacarrierIndexes-1, carriers(datacarrierIndexes));
stem(pilotTones-1, carriers(pilotTones),'g*');
hold off;
xlabel('Subcarrier Index');
ylabel('Magnitude');
xlim([0 nFFT-1]);
grid on;
legend('Guard Carriers','Data Carriers','Pilot Carriers');