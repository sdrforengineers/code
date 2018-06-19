%----------------------------------------------------------------------
% Chapter 2 
% "Digital Communication Systems Engineering Using Software Defined Radio
% MATLAB Scripts
%----------------------------------------------------------------------


% Clear workspace
clear all;close all;
addpath('support');

% Specify plot parameters
txtsize=10;
ltxtsize=9;
pwidth=4;
pheight=4;
pxoffset=0.65;
pyoffset=0.5;
markersize=5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
SetPlotSize ([pxoffset pyoffset pwidth pheight-1],'inches','white');
SetPlotFont ('Times', txtsize);
%
bandLimiting
%
set(gcf,'PaperPositionMode','auto');
print('bandLimiting','-depsc');
close;


%%
filteredEffectsBest
pause(1);
set(0, 'currentfigure', f1);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('filteredEffectsBest_a','-depsc');
pause(1);
close(f1);
set(0, 'currentfigure', f2);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('filteredEffectsBest_b','-depsc');
pause(1);
close(f2);
set(0, 'currentfigure', f3);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('filteredEffectsBest_c','-depsc');
pause(1);
close(f3);

%%
pause(1);
systemExample
pause(1);

set(0, 'currentfigure', f1);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('const_vs_time_a','-depsc');
close(f1);
pause(1);

set(0, 'currentfigure', f2);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('const_vs_time_b','-depsc');
close(f2);
pause(1);

set(0, 'currentfigure', f3);  %# for figures
SetPlotSize ([pxoffset pyoffset pwidth pheight-1],'inches','white');
SetPlotFont ('Times', txtsize);
set(gcf,'PaperPositionMode','auto');
print('const_vs_time_c','-depsc');
close(f3);
