x = logical([1 0 1 1 0 1 0 1 1 1 0 1]');
crcGen = comm.CRCGenerator('z^3 + 1');
crcDet = comm.CRCDetector('z^3 + 1');
codeword = crcGen(x);
codewordWithError = codeword; codewordWithError(1) = ~codewordWithError(1);
[tx, err]   = crcDet(codeword);
[tx1, err1] = crcDet(codewordWithError);