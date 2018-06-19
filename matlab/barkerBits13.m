% Show Barker Autocorrelations search example
sequenceLength = 13;
hBCode = comm.BarkerCode('Length',7,'SamplesPerFrame', sequenceLength);
seq = hBCode(); gapLen = 100; gapLenEnd = 200;
gen = @(Len) 2*randi([0 1],Len,1)-1;
y = [gen(gapLen); seq; gen(gapLenEnd)];
corr = xcorr(y,seq);
L = length(corr);
[v,i] = max(corr);
% Estimation of peak position
% The correlation sequence should be 2*L-1, where L is the length of the
% longest of the two sequences
%
% The first N-M will be zeros, where N is the length of the long sequence
% and N is the length of the shorter sequence
%
% The peak itself will occur at zero lag, or when they are directly
% overlapping, resulting in a peak in the middle sample in correlation
numZeros = length(y) - sequenceLength;
fromSeqEdge = gapLen + sequenceLength;
peakPosEst = numZeros + fromSeqEdge;
% Estimate 'gapLen' with known BarkerLength and Number of starting zeros
[~,p] = max(corr);
numZeros = length(y) - sequenceLength;
gapLenEst = p - numZeros - sequenceLength;