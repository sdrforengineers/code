%----------------------------------------------------------------------
% Chapter 6
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
barkerBits;
pause(2);
set(0, 'currentfigure', f1);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('barkerBits_a','-depsc');
close(f1);
pause(3);

set(0, 'currentfigure', f2);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('barkerBits_b','-depsc');
close(f2);
pause(2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compareBarkers;
pause(2);
set(0, 'currentfigure', f(1));  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('barkerAuto_a','-depsc');
close(f(1));
pause(3);
%
set(0, 'currentfigure', f(2));  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('barkerAuto_b','-depsc');
close(f(2));
pause(2);
%
set(0, 'currentfigure', f(3));  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('barkerAuto_c','-depsc');
close(f(3));
pause(2);
%
set(0, 'currentfigure', f(4));  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('barkerAuto_d','-depsc');
close(f(4));
pause(2);