%----------------------------------------------------------------------
% Chapter 5
% "Digital Communication Systems Engineering Using Software Defined Radio
% MATLAB Scripts"
%----------------------------------------------------------------------


% Clear workspace
clear all;
addpath('support');

% Specify plot parameters
txtsize=10;
ltxtsize=9;
pwidth=4;
pheight=4;
pxoffset=0.65;
pyoffset=0.5;
markersize=5;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqEstimate;
set(0, 'currentfigure', f1);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('sinusoids','-depsc');
close(f1);

set(0, 'currentfigure', f2);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('freqEst','-depsc');
close(f2);
pause(2);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqShiftFFT;
set(0, 'currentfigure', f1);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('PSD_nooffset','-depsc');
close(f1);

set(0, 'currentfigure', f2);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('PSD_offset','-depsc');
close(f2);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqCoarseEst;
pause(1);
set(0, 'currentfigure', f1);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('PSD2_nooffset','-depsc');
close(f1);
pause(1);

set(0, 'currentfigure', f2);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('PSD2_offset','-depsc');
close(f2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
viewRotation;
pause(1);
SetPlotSize ([pxoffset pyoffset pwidth-1 pheight-1],'inches','white');
SetPlotFont ('Times', txtsize+2);
set(gcf,'PaperPositionMode','auto');
print('rotatingConst','-depsc');
close;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plutoFreqOffset
pause(1);
set(0, 'currentfigure', f1);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth-1 pheight-1],'inches','white');
SetPlotFont ('Times', txtsize+2);
set(gcf,'PaperPositionMode','auto');
print('beforeTiming','-depsc');
close(f1);
pause(1);

set(0, 'currentfigure', f2);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth-1 pheight-1],'inches','white');
SetPlotFont ('Times', txtsize+2);
set(gcf,'PaperPositionMode','auto');
print('afterTiming','-depsc');
close(f2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badlock
pause(1);
%f1 = gca;
%set(0, 'currentfigure', f1);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth-1 pheight-1],'inches','white');
SetPlotFont ('Times', txtsize+2);
set(gcf,'PaperPositionMode','auto');
print('goodlock','-depsc');
close;
pause(1);

%f2 = gca;
%set(0, 'currentfigure', f2);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth-1 pheight-1],'inches','white');
SetPlotFont ('Times', txtsize+2);
set(gcf,'PaperPositionMode','auto');
print('badlock','-depsc');
close;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convergence;
pause(1);
set(0, 'currentfigure', f1);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('bloopLow','-depsc');
close(f1);
pause(1);

set(0, 'currentfigure', f2);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('bloopHigh','-depsc');
close(f2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estError;
pause(1);
SetPlotSize ([pxoffset pyoffset pwidth-0 pheight-1],'inches','white');
SetPlotFont ('Times', txtsize+2);
set(gcf,'PaperPositionMode','auto');
print('errorOverTime','-depsc');
close;



