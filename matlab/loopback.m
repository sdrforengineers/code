% Setup Receiver
rx=sdrrx('Pluto','OutputDataType','double','SamplesPerFrame',2^15);
% Setup Transmitter
tx = sdrtx('Pluto','Gain',-30);
% Transmit sinewave
sine = dsp.SineWave('Frequency',300,...
                    'SampleRate',rx.BasebandSampleRate,...
                    'SamplesPerFrame', 2^12,...
                    'ComplexOutput', true);
tx.transmitRepeat(sine()); % Transmit continuously
% Setup Scope
samplesPerStep = rx.SamplesPerFrame/rx.BasebandSampleRate;
steps = 3;
ts = dsp.TimeScope('SampleRate', rx.BasebandSampleRate,...
                   'TimeSpan', samplesPerStep*steps,...
                   'BufferLength', rx.SamplesPerFrame*steps);
% Receive and view sine
for k=1:steps
  ts(rx());
end
