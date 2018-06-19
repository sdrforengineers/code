for n = K+1:FrameLength
    % Apply equalizer to received data
    x_hat = f'*r(n:-1:n-K+1);
    % Estimate error
    if n<(PreambleLength-delta)
      e = x(i-delta) - x_hat;
    else
      e = sign(y) - x_hat;
    end
    % Update equalizer
    f = f + mu*conj(e)*r(n:-1:n-K+1);
end
