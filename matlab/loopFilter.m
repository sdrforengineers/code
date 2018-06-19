% Loop filter
loopFiltOut = LoopPreviousInput + LoopFilterState;
g = e*ProportionalGain + loopFiltOut; % Filter error signal
LoopFilterState = loopFiltOut;
LoopPreviousInput = e*IntegratorGain;
% Loop filter (alternative with filter objects)
lf = dsp.BiquadFilter('SOSMatrix',tf2sos([1 0],[1 -1])); % Create filter
g = lf(IntegratorGain*e) + ProportionalGain*e; % Filter error signal
