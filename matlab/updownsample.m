clc; clear all; clear time_plot; clear fft_plot; clear functions; file = 1;

% Specify plot parameters
txtsize=20;
ltxtsize=9;
pwidth=4; %4;
pheight=4; %4;
pxoffset=1; %0.65;
pyoffset=1; %0.5;
markersize=5;

% number filter taps
taps = 128;
% where to plot from
start = int32((taps/100)+1)*100;
inc = 50;

Fs1 = 10*2*pi;
% Create deterministic and stochastic digital data streams
n = 0:1/Fs1:100-(1/Fs1);                % Time index vector
sin_wave = sin(5*n*2*pi);               % Generation of sinusoidal signal
random = 2*round(rand(1,length(n)))-1;  % Random string of +1 and -1 values

figure(1);
plot(n(start:start+inc)*Fs1, sin_wave(start:start+inc), 'o', ...
    n(start:start+inc)*Fs1, sin_wave(start:start+inc));
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {'Original', 'Sinusoidal Signal : Time domain'});
figure(2)
plot(n(start:start+inc)*Fs1, random(start:start+inc),'o', ...
    n(start:start+inc)*Fs1, random(start:start+inc));
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {'Original', 'Random Binary Signal : Time domain'});

figure(3);
fft_plot(sin_wave(start:end), 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset,  pyoffset, markersize, ...
    {'Original', 'Sinusoidal Signal : Fourier domain'}, '$f_{s1}$');
figure(4);
fft_plot(random(start:end), 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {'Original', 'Random Binary Signal : Fourier domain'}, '$f_{s1}$');

% Create lowpass filter and apply it to both data streams
% b = firls(n,f,a),
%     n is the FIR filter order
%     f is a vector of pairs of frequency points,
%     a is a vector containing the desired amplitude at the points in f
coeffs1 = firls(taps,[0 0.2 0.22 1],[1 1 0 0]);   % FIR filter coefficients
sin_bwlimited = filter(coeffs1,1,sin_wave);
random_bwlimited = filter(coeffs1,1,random);

figure(5)
plot(n(start:start+inc)*Fs1, sin_bwlimited(start:start+inc), '-o');
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {'Band Limited', 'Sinusoidal Signal : Time domain'});
figure(6)
plot(n(start:start+inc)*Fs1, random_bwlimited(start:start+inc),'o');
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {'Band Limited', 'Random Binary Signal : Time domain'});

figure(7);
fft_plot(sin_bwlimited(start:end), 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {'Band Limited', 'Sinusoidal Signal : Fourier domain'}, '$f_{s1}$');
figure(8);
fft_plot(random_bwlimited(start:end), 1, txtsize, ltxtsize, pwidth, ...
    pheight, pxoffset, pyoffset, markersize, ...
    {'Band Limited', 'Random Binary Signal : Fourier domain'}, '$f_{s1}$');

% y = upsample(x,n)
%     increases the sampling rate of x by inserting (n – 1) zeros
%     between samples.
N = 5;
sin_up = upsample(sin_bwlimited,N);
random_up = upsample(random_bwlimited,N);
Fs2 = Fs1 * N;

figure(9)
plot(n(start*N:(start+inc)*N)*Fs1, sin_up(start*N:(start+inc)*N), 'o');
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d)', N), 'Sinusoidal Signal : Time domain'});
figure(10)
plot(n(start*N:(start+inc)*N)*Fs1, random_up(start*N:(start+inc)*N),'o');
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d)', N),  'Random Binary Signal : Time domain'});

figure(11);
fft_plot(sin_up(start:end), N, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d)', N), 'Sinusoidal Signal : Fourier domain'}, '$f_{s2}$');
figure(12);
fft_plot(random_up(start:end), N, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d)', N), 'Random Binary Signal : Fourier domain'}, '$f_{s2}$');

% Attempt to downsampling by M without filtering
% This is incorrect, but is instructive to show what artifacts occur
M = 3;
sin_up_down = downsample(sin_up,M);
random_up_down = downsample(random_up,M);
Fs3 = Fs2/M;

figure(13)
plot(n(start*N/M:(start+inc)*N/M)*Fs1, sin_up_down(start*N/M:(start+inc)*N/M), 'o');
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Improperly Upsampled (%d) then Downsampled (%d)', N, M), 'Sinusoidal Signal : Time domain'});
figure(14)
plot(n(start*N/M:(start+inc)*N/M)*Fs1, random_up_down(start*N/M:(start+inc)*N/M),'o');
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Improperly Upsampled (%d) then Downsampled (%d)', N, M), 'Random Binary Signal : Time domain'});

figure(15);
fft_plot(sin_up_down(start:end), N, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d) then Downsampled (%d)', N, M),  'Sinusoidal Signal : Fourier domain'}, '$f_{s3}$');
figure(16);
fft_plot(random_up_down(start:end), N, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d) then Downsampled (%d)', N, M), 'Random Binary Signal : Fourier domain'}, '$f_{s3}$');

% Lowpass filtering of baseband periodic replica followed by downsampling
% (correct approach)
coeffs2 = firls(taps,[0 0.15 0.25 1],[N N 0 0]); % FIR filter coefficients
sin_up_filtered = filter(coeffs2,1,sin_up);
sin_up_filtered_down = downsample(sin_up_filtered,M);
random_up_filtered = filter(coeffs2,1,random_up);
random_up_filtered_down = downsample(random_up_filtered,M);
start = start + int32(taps/(2*N));
figure(17)
plot(n(start*N:(start+inc)*N)*Fs1, sin_up_filtered(start*N:(start+inc)*N), 'o');
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Properly Upsampled (%d) and Filtered', N), 'Sinusoidal Signal : Time domain'});
figure(18)
plot(n(start*N:(start+inc)*N)*Fs1, random_up_filtered(start*N:(start+inc)*N),'o');
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Properly Upsampled (%d) and Filtered', N), 'Random Binary Signal : Time domain'});

figure(19);
fft_plot(sin_up_filtered(start:end), 1, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {sprintf('Properly Upsampled (%d) and Filtered', N),  'Sinusoidal Signal : Fourier domain'}, '$f_{s2}$');
figure(20);
fft_plot(random_up_filtered(start:end), 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Properly Upsampled (%d) and Filtered', N), 'Random Binary Signal : Fourier domain'}, '$f_{s2}$');

figure(21)
plot(n(start*N/M:(start+inc)*N/M)*Fs1, sin_up_filtered_down(start*N/M:(start+inc)*N/M), 'o');
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Properly Upsampled (%d), Filtered, and Downsampled (%d)', N, M), 'Sinusoidal Signal : Time domain'});
figure(22)
plot(n(start*N/M:(start+inc)*N/M)*Fs1, random_up_filtered_down(start*N/M:(start+inc)*N/M),'o');
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Properly Upsampled (%d), Filtered, and Downsampled (%d)', N, M), 'Random Binary Signal : Time domain'});

figure(23);
fft_plot(sin_up_filtered_down(start:end), 1, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {sprintf('Properly Upsampled (%d), Filtered and Downsampled (%d)', N, M),  'Sinusoidal Signal : Fourier domain'}, '$f_{s3}$');
figure(24);
fft_plot(random_up_filtered_down(start:end), 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Properly Upsampled (%d), Filtered and Downsampled (%d)', N, M), 'Random Binary Signal : Fourier domain'}, '$f_{s3}$');

clear time_plot; clear fft_plot; clear functions;
Vars=whos;
PersistentVars=Vars([Vars.persistent]);
PersistentVarNames={PersistentVars.name};
clear(PersistentVarNames{:});

return;
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
function fft_plot(data, points, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, titlestr, units)
    persistent file;
    fft_data = do_fft(data);
    fft_axis=do_fs(data);
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
    xlim([-0.5 0.5]);
    xticks([-.5 -.25  0  .25  .5]);
    grid on;
    title(titlestr);
end
