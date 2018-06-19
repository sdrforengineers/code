%----------------------------------------------------------------------
% Chapter 3 
% "Digital Communication Systems Engineering Using Software Defined Radio
% MATLAB Scripts
%----------------------------------------------------------------------

% Clear workspace
clear all;


% Specify plot parameters
txtsize=10;
ltxtsize=9;
pwidth=4;
pheight=4;
pxoffset=0.65;
pyoffset=0.5;
markersize=5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating a binary data source

% Create three binary data streams of length L, each with a different percentage of 1
% and 0 values
L = 100; 
prob1 = 0.5; prob2 = 0.6; prob3 = 0.2; % Probability values for 1 outputs
b1 = round(0.5*rand(1,L)+0.5*prob1); % Even split between 1 and 0 values
b2 = round(0.5*rand(1,L)+0.5*prob2); % Have 60% 1 values and 40% 0 values
b3 = round(0.5*rand(1,L)+0.5*prob3); % Have 20% 1 values and 80% 0 values

% Display data stream and associated histogram
figure(1);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
stem((0:1:L-1),b1);xlabel('Discrete Time (n)');ylabel('Digital Value');
%title(sprintf('Binary Data Stream\n (%d%% ones, %d%% zeros)',100*(prob1),100*(1-prob1)));
set(gcf,'PaperPositionMode','auto');
print('../ch3_binprob_s1','-depsc');
close;
figure(1);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
histogram(b1,2);xlabel('Digital Value');ylabel({'Occurance';'Probability (%)'});
%title(sprintf('Histogram of Binary Stream\n (%d%% ones, %d%% zeros)',100*(prob1),100*(1-prob1)));
set(gcf,'PaperPositionMode','auto');
print('../ch3_binprob_h1','-depsc');
close;
figure(1);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
stem((0:1:L-1),b2);xlabel('Discrete Time (n)');ylabel('Digital Value');
%title(sprintf('Binary Data Stream\n (%d%% ones, %d%% zeros)',100*(prob2),100*(1-prob2)));
set(gcf,'PaperPositionMode','auto');
print('../ch3_binprob_s2','-depsc');
close;
figure(1);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
histogram(b2,2);xlabel('Digital Value');ylabel({'Occurance';'Probability (%)'});
%title(sprintf('Histogram of Binary Stream\n (%d%% ones, %d%% zeros)',100*(prob2),100*(1-prob2)));
set(gcf,'PaperPositionMode','auto');
print('../ch3_binprob_h2','-depsc');
close;
figure(1);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
stem((0:1:L-1),b3);xlabel('Discrete Time (n)');ylabel('Digital Value');
%title(sprintf('Binary Data Stream\n (%d%% ones, %d%% zeros)',100*(prob3),100*(1-prob3)));
set(gcf,'PaperPositionMode','auto');
print('../ch3_binprob_s3','-depsc');
close;
figure(1);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
histogram(b3,2);xlabel('Digital Value');ylabel({'Occurance';'Probability (%)'});
%title(sprintf('Histogram of Binary Stream\n (%d%% ones, %d%% zeros)',100*(prob3),100*(1-prob3)));
set(gcf,'PaperPositionMode','auto');
print('../ch3_binprob_h3','-depsc');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binary channel example

% Define simulation parameters
L = 100000; % Transmission length
prob00 = 0.95; % Probability zero received given zero transmitted
prob11 = 0.99; % Probability one received given one transmitted
prob4 = 0.7; % Have 40% 1 values and 60% 0 values

% Create transmitted binary data stream
b4 = round(0.5*rand(1,L)+0.5*prob4); % Have 40% 1 values and 60% 0 values
b4hat = b4; % Initialize receive binary data stream

% Randomly select 1 and 0 values for flipping
ind_zero = find(b4 == 0); % Find those time instances with zero values
ind_one = find(b4 == 1); % Find those time instances with one values
ind_flip_zero = find(round(0.5*rand(1,length(ind_zero))+0.5*(1-prob00)) == 1); % Flip zero bits to one bits
ind_flip_one = find(round(0.5*rand(1,length(ind_one))+0.5*(1-prob11)) == 1); % Flip one bits to zero bits

% Corrupt received binary data stream
b4hat(ind_zero(ind_flip_zero)) = 1; % Flip 0 to 1
b4hat(ind_one(ind_flip_one)) = 0; % Flip 1 to 0

% Calculate bit error statistics
b4error_total = sum(abs(b4-b4hat))/L;
b4error_1 = sum(abs(b4(ind_one) - b4hat(ind_one)))/length(ind_one);
b4error_0 = sum(abs(b4(ind_zero) - b4hat(ind_zero)))/length(ind_zero);

% Plot histogram of bit errors
figure(2);
SetPlotSize ([pxoffset pyoffset pwidth pheight-2.5],'inches','white');
SetPlotFont ('Times', txtsize);
bar([1 2 3],[b4error_total b4error_1 b4error_0]);xlabel('Bit Error Rates: (1) Total, (2) One Transmitted, (3) Zero Transmitted');ylabel('Probability');
%title('Binary Channel Error Probabilities');
set(gcf,'PaperPositionMode','auto');
print('../ch3_binchannelerr','-depsc');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Random Variable PDFs and CDFs

% Define simulation parameters
L = 1000000; % Length of random samples
mean_normal = 1; stddev_normal = 2; % Mean and standard deviation of normal RV values
res_hist = 100; % Histogram resolution

% Generate random samples for different distributions
b_unif = rand(1,L); % Uniform random values
b_normal = (stddev_normal*randn(1,L) + mean_normal); % Normal random values
b_rayleigh = (sqrt(randn(1,L).^2 + randn(1,L).^2)); % Rayleigh random values

% Obtain cumulative distribution functions
[N_unif, edges_unif] = histcounts(b_unif,res_hist);
N_unif_cum = cumsum(N_unif)./L;
[N_normal, edges_normal] = histcounts(b_normal,res_hist);
N_normal_cum = cumsum(N_normal)./L;
[N_rayl, edges_rayl] = histcounts(b_rayleigh,res_hist);
N_rayl_cum = cumsum(N_rayl)./L;

% Calculate probability of values between 0.7 and 1
x_lower = 0.7;
x_upper = 1.0;
unif_ind_range = find((x_lower <= edges_unif) & (edges_unif < x_upper));
normal_ind_range = find((x_lower <= edges_normal) & (edges_normal < x_upper));
rayl_ind_range = find((x_lower <= edges_rayl) & (edges_rayl < x_upper));
prob_unif = sum(N_unif(unif_ind_range))./L;
prob_normal = sum(N_normal(normal_ind_range))./L;
prob_rayl = sum(N_rayl(rayl_ind_range))./L;

% Plot the PDF and CDF of these random values
figure(3);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
histogram(b_unif,res_hist,'Normalization','probability');xlabel('Output Values of Random Variable');ylabel('Probability');
%title({'Uniform Probability';'Density Function'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_pdfs_cdfs_h1','-depsc');
close;
figure(3);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
bar((edges_unif(1:1:(length(edges_unif)-1))+edges_unif(2:1:(length(edges_unif))))./2,N_unif_cum,(edges_unif(2)-edges_unif(1)));xlabel('Output Values of Random Variable');ylabel({'Cumulative';'Probability'});hold on;line([x_lower x_lower],[0 1]);line([x_upper x_upper],[0 1]);text(0.4,0.9,sprintf('P(%0.5g<X<%0.5g)=%0.5g',x_lower,x_upper,prob_unif));hold off;
%title({'Uniform Cumulative';'Distribution Function'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_pdfs_cdfs_b1','-depsc');
close;
figure(3);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
histogram(b_normal,res_hist,'Normalization','probability');xlabel('Output Values of Random Variable');ylabel('Probability');
%title({'Gaussian Probability';'Density Function'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_pdfs_cdfs_h2','-depsc');
close;
figure(3);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
bar((edges_normal(1:1:(length(edges_normal)-1))+edges_normal(2:1:(length(edges_normal))))./2,N_normal_cum,(edges_normal(2)-edges_normal(1)));xlabel('Output Values of Random Variable');ylabel({'Cumulative';'Probability'});hold on;line([x_lower x_lower],[0 1]);line([x_upper x_upper],[0 1]);text(-5,0.9,sprintf('P(%0.5g<X<%0.5g)=%0.5g',x_lower,x_upper,prob_normal));hold off;
%title({'Gaussian Cumulative';'Distribution Function'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_pdfs_cdfs_b2','-depsc');
close;
figure(3);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
histogram(b_rayleigh,res_hist,'Normalization','probability');xlabel('Output Values of Random Variable');ylabel('Probability');
%title({'Rayleigh Probability';'Density Function'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_pdfs_cdfs_h3','-depsc');
close;
figure(3);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
bar((edges_rayl(1:1:(length(edges_rayl)-1))+edges_rayl(2:1:(length(edges_rayl))))./2,N_rayl_cum,(edges_rayl(2)-edges_rayl(1)));xlabel('Output Values of Random Variable');ylabel({'Cumulative';'Probability'});hold on;line([x_lower x_lower],[0 1]);line([x_upper x_upper],[0 1]);text(0.1,0.9,sprintf('P(%0.5g<X<%0.5g)=%0.5g',x_lower,x_upper,prob_rayl));hold off;
%title({'Rayleigh Cumulative';'Distribution Function'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_pdfs_cdfs_b3','-depsc');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate bivariate random variables 

% Define simulation parameters
L = 1000000; % Length of data streams
res_hist = 100; % Histogram resolution
std_dev = 5; % Standard deviation of input Gaussian variables

% Create uncorrelated 2D Gaussian random variable
x_normal_1 = std_dev.*randn(1,L);
y_normal_1 = std_dev.*randn(1,L);

% Create correlated 2D Gaussian random data stream
x_normal_2 = x_normal_1+0.1*y_normal_1;
y_normal_2 = y_normal_1+0.9*x_normal_1;

% Plot 2D histograms of uncorrelated and correlated 2D Gaussian random data stream
figure(4);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
histogram2(x_normal_1,y_normal_1,res_hist,'Normalization','pdf','DisplayStyle','bar3');xlabel('X');ylabel('Y');
%title('Uncorrelated Gauassian');
set(gcf,'PaperPositionMode','auto');
print('../ch3_bivargauss_uncorrb','-depsc');
close;
figure(4);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
histogram2(x_normal_1,y_normal_1,res_hist,'Normalization','pdf','DisplayStyle','tile');xlabel('X');ylabel('Y');
%title('Uncorrelated Gauassian');
set(gcf,'PaperPositionMode','auto');
print('../ch3_bivargauss_uncorrt','-depsc');
close;
figure(4);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
histogram2(x_normal_2,y_normal_2,res_hist,'Normalization','pdf','DisplayStyle','bar3');xlabel('X');ylabel('Y');
%title('Correlated Gaussian');
set(gcf,'PaperPositionMode','auto');
print('../ch3_bivargauss_corrb','-depsc');
close;
figure(4);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-2],'inches','white');
SetPlotFont ('Times', txtsize);
histogram2(x_normal_2,y_normal_2,res_hist,'Normalization','pdf','DisplayStyle','tile');xlabel('X');ylabel('Y');
%title('Correlated Gaussian');
set(gcf,'PaperPositionMode','auto');
print('../ch3_bivargauss_corrt','-depsc');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observe impact of filtering random processes 

% Define simulation parameters
L = 1000000; % Length of data streams
res_hist = 100; % Histogram resolution
cutoff_freq1 = 0.2; % Small passband LPF
cutoff_freq2 = 0.7; % large passband LPF
filt_coeffs1 = firls(13,[0 cutoff_freq1 cutoff_freq1+0.02 1],[1 1 0 0]);
filt_coeffs2 = firls(13,[0 cutoff_freq2 cutoff_freq2+0.02 1],[1 1 0 0]);

% Create input 2D Gaussian random variable
x_in = rand(1,L);
y_in = rand(1,L);

% Filter input random data stream
filt_output1 = filter(filt_coeffs1,1,(x_in+1j.*y_in));
x_out1 = real(filt_output1);
y_out1 = imag(filt_output1);
filt_output2 = filter(filt_coeffs2,1,(x_in+1j.*y_in));
x_out2 = real(filt_output2);
y_out2 = imag(filt_output2);

% Plot 2D histogram of input and output random data stream
figure(5);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
histogram2(x_in,y_in,res_hist,'Normalization','pdf','DisplayStyle','bar3');xlabel('X');ylabel('Y');
%title({'Probability';'Input Random Process'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_psdfilt_inb','-depsc');
close;
figure(5);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
histogram2(x_in,y_in,res_hist,'Normalization','pdf','DisplayStyle','tile');xlabel('X');ylabel('Y');
%title({'Probability';'Input Random Process'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_psdfilt_int','-depsc');
close;
figure(5);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
histogram2(x_out1,y_out1,res_hist,'Normalization','pdf','DisplayStyle','bar3');xlabel('X');ylabel('Y');
%title({'Probability';'Output Random Process 1'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_psdfilt_out1b','-depsc');
close;
figure(5);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
histogram2(x_out1,y_out1,res_hist,'Normalization','pdf','DisplayStyle','tile');xlabel('X');ylabel('Y');
%title({'Probability';'Output Random Process 1'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_psdfilt_out1t','-depsc');
close;
figure(5);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
histogram2(x_out2,y_out2,res_hist,'Normalization','pdf','DisplayStyle','bar3');xlabel('X');ylabel('Y');
%title({'Probability';'Output Random Process 2'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_psdfilt_out2b','-depsc');
close;
figure(5);
SetPlotSize ([pxoffset pyoffset pwidth-2 pheight-3],'inches','white');
SetPlotFont ('Times', txtsize);
histogram2(x_out2,y_out2,res_hist,'Normalization','pdf','DisplayStyle','tile');xlabel('X');ylabel('Y');
%title({'Probability';'Output Random Process 2'});
set(gcf,'PaperPositionMode','auto');
print('../ch3_psdfilt_out2t','-depsc');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% END OF CHAPTER 3 EXAMPLES %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
