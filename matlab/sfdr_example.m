clear all;

deltat = 1e-8;
fs = 1/deltat;
t = 0:deltat:1e-5-deltat;
fundamental = 3959297;  % Prime number
x = 10e-3*sin(2*pi*fundamental*t);
sfdr(x,fs);