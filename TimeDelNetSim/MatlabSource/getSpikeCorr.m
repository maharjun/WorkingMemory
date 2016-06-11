function [ CorrCoeffs ] = getSpikeCorr( SpikeMat1, SpikeMat2, Jitter)
%GETSPIKECORR Returns correlation of two Spike Matrices
%   Detailed explanation goes here

FullSM1 = full(SpikeMat1);
FullSM2 = full(SpikeMat2);

[FilterSL1, z1] = filter(triang(2*Jitter+1), 1, FullSM1, [], 2);
[FilterSL2, z2] = filter(triang(2*Jitter+1), 1, FullSM2, [], 2);

FilterSL1 = [FilterSL1(:, Jitter+1:end), z1(1:Jitter, :)'];
FilterSL2 = [FilterSL2(:, Jitter+1:end), z2(1:Jitter, :)'];

Length1 = size(FullSM1); Length1 = Length1(2);
Length2 = size(FullSM2); Length2 = Length2(2);

LengthTot = Length1+Length2-1;
Min2Pow   = 2^ceil(log2(LengthTot));

FFTFiltSL1 = fft(FilterSL1, Min2Pow, 2);
FFTFiltSL2 = fft(FilterSL2, Min2Pow, 2);

FFTCorr = FFTFiltSL2.*conj(FFTFiltSL1);

CorrMat = ifft(FFTCorr, [], 2, 'symmetric');
CorrMat = CorrMat(:, 1:LengthTot);

CorrCoeffs = sum(CorrMat);
end


