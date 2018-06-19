% ZC-TED calculation occurs on a strobe
if Trigger && all(~TriggerHistory(2:end))
    % Calculate the midsample point for odd or even samples per symbol
    t1 = TEDBuffer(end/2 + 1 - rem(N,2));
    t2 = TEDBuffer(end/2 + 1);
    midSample = (t1+t2)/2;
    e = real(midSample)*(sign(real(TEDBuffer(1)))-sign(real(filtOut))) ...
        imag(midSample)*(sign(imag(TEDBuffer(1)))-sign(imag(filtOut)));
else
    e = 0;
end
% Update TED buffer to manage symbol stuffs
switch sum([TriggerHistory(2:end), Trigger])
    case 0
      % No update required
    case 1
      % Shift TED buffer regularly if ONE trigger across N samples
      TEDBuffer = [TEDBuffer(2:end), filtOut];
    otherwise % > 1
      % Stuff a missing sample if TWO triggers across N samples
      TEDBuffer = [TEDBuffer(3:end), 0, filtOut];
end