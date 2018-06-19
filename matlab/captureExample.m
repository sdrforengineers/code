% Transmit frame repeatedly
tx = sdrtx('Pluto');
tx = sdrtx('Pluto','SamplesPerFrame',length(frame)*2);
tx.transmitRepeat(frame);
for k=1:4,rx();end; % Remove stale data from buffers
rxBuffer = rx();
