% Interpolation Controller with modulo-1 counter
d = g + 1/N;
TriggerHistory = [TriggerHistory(2:end), Trigger];
Trigger = (Counter < d); % Check if a trigger condition
if Trigger % Upate mu if a trigger
    mu = Counter / d;
end
Counter = mod(Counter - d, 1); % Update counter