% Define interpolator coefficients
alpha = 0.5;
InterpFilterCoeff = ...
    [ 0,       0,         1,       0;    % Constant
    -alpha, 1+alpha, -(1-alpha), -alpha; % Linear
    alpha,  -alpha,    -alpha,   alpha]; % Quadratic
% Filter input data
ySeq = [y(i); InterpFilterState]; % Update delay line
 % Produce filter output
filtOut = sum((InterpFilterCoeff * ySeq) .* [1; mu; mu^2]);
InterpFilterState = ySeq(1:3); % Save filter input data