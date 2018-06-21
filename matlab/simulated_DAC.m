clc; clear all; clear time_plot; clear fft_plot; clear functions;

% Specify plot parameters
txtsize=10;
ltxtsize=9;
pwidth=4;
pheight=4;
pxoffset=0.65;
pyoffset=0.5;
markersize=5;
inc=20;
start = 1;
bits=3;
mult = 10;

Fs1 = 10*2*pi;
% Create data streams
tanalog = 0:1/(Fs1*mult):100-(1/(Fs1*mult));   % Time index vector
tdigital = 0:1/(Fs1):100-(1/Fs1);
sin_analog = sin(5*tanalog*2*pi);         % Generation of sinusoidal signal
sin_digital = sin(5*tdigital*2*pi);
sin_quantized = (floor(sin_analog*(2^bits)))/(2^bits);
sin_sample = downsample(sin_quantized, mult);
sin_time_sample = zeros(1,length(tdigital));

for i = 0:length(sin_sample)-1;
    for j = 1:mult
        sin_time_sample(i*mult + j) = sin_sample(i+1);
    end
end

figure(1);
plot(tdigital(1:inc)*Fs1, sin_digital(1:inc), 'o', ...
    tanalog(1:inc*mult)*Fs1, sin_analog(1:inc*mult));
time_plot(1, inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {'Analog', 'Sinusoidal Signal : Time domain'});

figure(2);
plot(tanalog(1:inc*mult)*Fs1, sin_quantized(1:inc*mult), ...
    tdigital(1:inc)*Fs1, sin_sample(1:inc), 'o');
time_plot(1, inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {'Quanitized', 'Sinusoidal Signal : Time domain'});

figure(3);
stairs(tdigital(1:inc)*Fs1, sin_sample(1:inc), '-o');
time_plot(1, inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {'Sampled', 'Sinusoidal Signal : Time domain'});

figure(4);
fft_plot(sin_time_sample, 1, mult, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset,  pyoffset, markersize, ...
    {'Sampled', 'Sinusoidal Signal : Fourier domain'}, '$f_s$');

% takes in of time domain vector, and provides
% output of FFT in dBFS (full scale)
function out = do_fft(in)
    L = length (in); % Window Length of FFT
    in_HannWnd = in' .* hanning(L ,'periodic');
    out = 20*(log10(abs(fftshift(fft(in_HannWnd,L))/(L/2))));
end
function out = do_fs(in)
   out = ((0:1:(length(in)-1))/(0.5*length(in))-1)/2;
end
function time_plot(x1, x2, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, titlestr)
    persistent file;
    xlabel('Discrete Time (n)');ylabel('Signal Amplitude');
    ylim([-1.5 1.5]);
    xlim([x1 x2]);
    title(titlestr);
end
function fft_plot(data, points, oversample, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, titlestr, units)
    persistent file;
    fft_data = do_fft(data);
    fft_axis=do_fs(data)*oversample;
    [psor,lsor] = findpeaks(fft_data(int32(length(data)/2):length(data)), ...
        fft_axis(int32(length(data)/2):length(data)), ...
    'NPeaks', points, 'SortStr', 'descend');
    if mean(psor) > -23
        [M,I] = min(lsor);
        plot(fft_axis, fft_data, lsor(I), psor(I), 'o');
        xlabel('$\frac{f_s}{2}$','interpreter','latex');
        str = sprintf('%2.1fdB @ %.3f', psor(I), lsor(I));
        text(lsor(I)+.02,psor(I), ...
            strcat(str, '$\frac{f_s}{2}$'), ...
            'interpreter','latex');
    else
        plot(fft_axis, fft_data);
    end
    xlabel(strcat({'Frequency'}, {' '}, {units}), 'interpreter','latex');
    ylabel('Signal Amplitude (dB)');
    ylim([-200 0]);
    xlim([-0.5*oversample 0.5*oversample]);
    xticks(-oversample/2:.5:oversample/2);
    grid on;
    title(titlestr);
end
