function [S Cx] = cyclic(x)
x = x';
y = x;
%N must be even and divisible by 4 and < lx
N = 512;
lx=length(x);
%These should be tweaked during network start up
%for best performance.  The worse the channel, the lower they should be.
%In an ideal channel, they should be closer to 0.8

%Set up variables
n=0:floor(lx-N);
ln=length(n);
%Compute windowing functions for later.
a=feval('hamming',N)';
g=feval('hamming',ln)';
g=g/sum(g);
a=a/sum(a);
Ts=1/N;
%Pre-allocate for speed
S=zeros(N+1,N/2+1);
X=zeros(2*N+1,ln);
Y=zeros(2*N+1,ln);

%Freq. Smoothed Cyclic Periodogram
for f=-N:N
  %N point FFTs of the signal are computed
  xf=x.*exp(-1i*2*pi*f*(0:lx-1)*Ts);
  yf=y.*exp(-1i*2*pi*f*(0:lx-1)*Ts);
  for i=1:ln
    %Multiply the FFT of X with the conj of Y and vice versa
    n_r=n(i)+(1:N);
    X(f+N+1,i)=a*xf(n_r)';
    Y(f+N+1,i)=conj(a*yf(n_r)');	
  end
end

for alpha=-N/4:N/4
  for f=-N/2:N/2
    f1=f+alpha;
    f2=f-alpha;
    if (abs(f1)<N/2)&&(abs(f2)<N/2)
      %g acts to smooth X*Y out, this is more obvious if you plot g
      %s is the cross correlation of X's and Y's frequency components
      %seperated by f +/- alpha
      S(f+N/2+1,N/4+alpha+1)=g*(X(f1+N+1,:).*Y(f2+N+1,:))';
    end
  end
end

%Compute correlation coefficients
Cx = fftshift(corrcoef(S').^2);
end