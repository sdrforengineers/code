%Estimate BER for Rate 2/3 Convolutional Code
%This example performs a bit error rate simulation for a link that uses 16-QAM modulation a rate 2/3 convolutional code.
%Set the modulation order, and compute the number of bits per symbol.
M = 16;
k = log2(M);
%Define a convolutinal coding trellis for a rate 2/3 code.
tPoly = poly2trellis([5 4],[23 35 0; 0 5 13]);
codeRate = 2/3;


EbNos = 0:0.3:1;
bers = zeros(40,length(EbNos));

for tblt=5:40
    for Index = 1:length(EbNos)
        berm = zeros(1e2,1);
        for g=1:100
            
            %Generate random binary data.
            dataIn = randi([0 1],100000,1);
            
            
            %Convolutinally encode the input data.
            codeword = convenc(dataIn,tPoly);
            
            %Reshape the encoded column vector into a matrix having k columns. Then, convert the binary matrix into an integer column vector.
            codewordMat = reshape(codeword,length(codeword)/k,k);
            txSym = bi2de(codewordMat);
            
            %Apply 16-QAM modulation to the encoded symbols.
            txSig = qammod(txSym,M);
            
            %Convert a 10 dB Eb/No to an equivalent signal-to-noise ratio. Pass the signal through an AWGN channel.
            EbNo = EbNos(Index);
            snr = EbNo + 10*log10(k*codeRate);
            rxSig = awgn(txSig,snr,'measured');
            
            demodSig = qamdemod(rxSig,M);
            
            demodSigMat = de2bi(demodSig,k);
            demodSigBinary = demodSigMat(:);
            
            traceBack = 16;
            
            dataOut = vitdec(demodSigBinary,tPoly,tblt,'cont','hard');
            
            decDelay = 2*traceBack;
            [numErrors,ber] = biterr(dataIn(1:end-decDelay),dataOut(decDelay+1:end));
            berm(g) = ber;
            %berUncoded = berawgn(EbNo,'qam',M);
            %berUncoded/ber
            %The convolutional code reduces the BER by approximately a factor of 4.
        end
        disp(['EbNo: ',num2str(EbNos(Index))]);
        bers(tblt,Index) = mean(berm);
    end
end

% Plot
bers = bers(5:40,:);
bersm = mean(bers,2);
plot(EbNos,bersm);xlabel('EbN0');ylabel('TBL');


