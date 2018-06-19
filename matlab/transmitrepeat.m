% Transmit all zeros
tx = sdrtx('Pluto');
fs = 1e6; fc = 1e4; s = 2*pi*fs*fc*(1:2^14).';
wave = complex(cos(s),sin(s));
tx.transmitRepeat(wave);
