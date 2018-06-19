% View some spectrum
rx = sdrrx('Pluto');
rx.SamplesPerFrame = 2^15;
sa = dsp.SpectrumAnalyzer;
sa.SampleRate = rx.BasebandSampleRate;
for k=1:1e3
  sa(rx());
end
