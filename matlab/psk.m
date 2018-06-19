clc; clear all; clear time_plot; clear fft_plot; clear functions;

% Specify plot parameters
txtsize=10;
ltxtsize=9;
pwidth=4;
pheight=4;
pxoffset=0.65;
pyoffset=0.5;
markersize=5;
fs=2*pi*100;
t=0:1/fs:100-(1/fs);
f1=20; %input(‘Carrier Sine wave frequency =’);
f2=5; %input(‘Message frequency =’);
x=sin(f1*t);%Carrier Sine

figure(1);
plot(t*fs,x);

time_plot(1000, 3000, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    'Carrier : Time domain');

figure(2);
fft_plot(x, 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset,  pyoffset, markersize, ...
    {'Carrier : Fourier domain'}, '$f_s$');

%Message signal
u=square(f2*t);

figure(3);
plot(t*fs,u);
time_plot(1000, 3000,  txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    'Binary Message');

figure(4);
fft_plot(u, 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset,  pyoffset, markersize, ...
    {'Binary Message : Fourier domain'}, '$f_s$');


v=x.*u;%Sine wave multiplied with square wave

figure(5);
plot(t*fs,v);
time_plot(1000, 3000, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    'PSK Modulated Waveform');

figure(6);
fft_plot(u, 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset,  pyoffset, markersize, ...
    {'PSK Modulated Waveform', ' domain'}, '$f_{s1}$');


return
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
    save=1;
    if save
        if isempty(file)  || (file >= 13)
            file = 1;
        else
            file = file + 1;
        end
        SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
        SetPlotFont ('Times', txtsize);
        set(gcf,'PaperPositionMode','auto');
        print(sprintf('ch2_psk_time_%d', file),'-depsc');
    end
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
    ylim([-100 0]);
    xlim([-0.5 0.5]);
    xticks([]);
    grid on;
    save=1;
    if save
        if isempty(file) ||  (file >= 13)
            file = 1;
        else
            file = file + 1;
        end
        SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
        SetPlotFont ('Times', txtsize);
        set(gcf,'PaperPositionMode','auto');
        print(sprintf('ch2_psk_fft_%d', file),'-depsc');
    end
    title(titlestr);
end