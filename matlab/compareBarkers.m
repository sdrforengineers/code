% Plot Barker sequences
%C = [1, 2, 3, 4, 5, 7, 11, 13];
C = [5, 7, 11, 13]; f = zeros(size(C));
for k=1:length(C)
    hBCode = comm.BarkerCode('Length',C(k),'SamplesPerFrame', C(k));
    seq = step(hBCode);
    f(k) = figure;
    %subplot(length(C),1,k);
    cc = xcorr(seq);
    stem(cc);
    grid on;ylim([-1 21]);xlim([0 C(end)*2-1]);
    ylabel('Autocorr');xlabel('Samples');
    %title(['N=',num2str(C(k))]);
    [~,i] = max(xcorr(seq));
    disp([num2str(C(k)),' ',num2str(i)]);
end