% Specify plot parameters
txtsize=10;
ltxtsize=9;
pwidth=4;
pheight=4;
pxoffset=0.65;
pyoffset=0.5;
markersize=5;

deltat = 1e-8;
fs = 1/deltat;
t = 0:deltat:1e-5-deltat;
fundamental=3959297;
x = 10e-3*sin(2*pi*fundamental*t);
r = sfdr(x,fs)
f1=figure(1);
sfdr(x,fs);
ylim([-400 10]);
f2=figure(2);
plot(t*10e6,x);
xlabel('time (us)');ylabel('Amplitude');
ylim([-1.1 1.1]);

%% GENERATE FIGURES HERE
figure(1);
set(0, 'currentfigure', f1);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_float.eps','-depsc');

figure(2);
set(0, 'currentfigure', f2);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_float_time.eps','-depsc');

bits=2^11;
x = round(10e-3*bits*sin(2*pi*fundamental*t))/bits;
r = sfdr(x,fs)
f1=figure(1);
sfdr(x,fs);
ylim([-400 10]);
f2=figure(2);
plot(t*10e6,x);
xlabel('time (us)');ylabel('Amplitude');
ylim([-1.1 1.1]);

%% GENERATE FIGURES HERE
figure(1);
set(0, 'currentfigure', f1);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_small.eps','-depsc');

figure(2);
set(0, 'currentfigure', f2);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_time_small.eps','-depsc');


bits=2^11;
x = round(bits*sin(2*pi*fundamental*t))/bits;
r = sfdr(x,fs)
f1=figure(1);
sfdr(x,fs);
ylim([-400 10]);
f2=figure(2);
plot(t*10e6,x);
xlabel('time (us)');ylabel('Amplitude');
ylim([-1.1 1.1]);

%% GENERATE FIGURES HERE
figure(1);
set(0, 'currentfigure', f1);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit.eps','-depsc');

figure(2);
set(0, 'currentfigure', f2);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_time.eps','-depsc');

t = 0:deltat:1e-4-deltat;
x = round(bits*sin(2*pi*fundamental*t))/bits;
r = sfdr(x,fs)
f1=figure(1);
sfdr(x,fs);
ylim([-400 10]);
f2=figure(2);
plot(t*10e6,x);
xlabel('time (us)');ylabel('Amplitude');
ylim([-1.1 1.1]);

%% GENERATE FIGURES HERE
figure(1);
set(0, 'currentfigure', f1);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_l.eps','-depsc');

figure(2);
set(0, 'currentfigure', f2);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_time_l.eps','-depsc');

fundamental=4000000;
x = round(bits*sin(2*pi*fundamental*t))/bits;
r = sfdr(x,fs)
f1=figure(1);
sfdr(x,fs);
ylim([-400 10]);
f2=figure(2);
plot(t*10e6,x);
xlabel('time (us)');ylabel('Amplitude');
ylim([-1.1 1.1]);

%% GENERATE FIGURES HERE
figure(1);
set(0, 'currentfigure', f1);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_corr.eps','-depsc');

figure(2);
set(0, 'currentfigure', f2);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_time_corr.eps','-depsc');

ran = rand(1,length(t)) - 0.5;
x = round(bits*sin(2*pi*fundamental*t) + ran)/bits;
r = sfdr(x,fs)
f1=figure(1);
sfdr(x,fs);
ylim([-400 10]);
f2=figure(2);
plot(t*10e6,x);
xlabel('time (us)');ylabel('Amplitude');
ylim([-1.1 1.1]);

%% GENERATE FIGURES HERE
figure(1);
set(0, 'currentfigure', f1);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_dither.eps','-depsc');

figure(2);
set(0, 'currentfigure', f2);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_time_dither.eps','-depsc');

fundamental=3959297;
x = round(bits*sin(2*pi*fundamental*t) + ran)/bits;
r = sfdr(x,fs)
f1=figure(1);
sfdr(x,fs);
ylim([-400 10]);
f2=figure(2);
plot(t*10e6,x);
xlabel('time (us)');ylabel('Amplitude');
ylim([-1.1 1.1]);

%% GENERATE FIGURES HERE
figure(1);
set(0, 'currentfigure', f1);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_dither_corr.eps','-depsc');

figure(2);
set(0, 'currentfigure', f2);  % Optional select given figure from handle
%%%
SetPlotSize ([pxoffset pyoffset pwidth pheight],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('./ch2_sfdr_12bit_time_dither_corr.eps','-depsc');