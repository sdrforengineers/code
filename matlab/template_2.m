%% Template 2
% Load data and perform processing
bfr = comm.BasebandFileReader(bfw.Filename, 'SamplesPerFrame',frameSize);
sa2 = dsp.SpectrumAnalyzer;
% Process each frame from the saved file
for frame = 1:framesToCollect
    sa2(bfr()); % Algorithm processing
end
