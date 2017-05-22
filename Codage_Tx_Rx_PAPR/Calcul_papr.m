function [paprSymboldB] = Calcul_papr(input_OFDM, nSymbol)

deb = 1; fin = nSymbol;
nBitperSymbol = 80;  % nFFTSize(64) + n_IG(16) = 80

input_OFDM = reshape(input_OFDM, nBitperSymbol, nSymbol); 

for ii = deb:fin;
    
    x = transpose(input_OFDM(:,ii));
        
    % Calcul du peak to average power ratio PAPR
    meanSquareValue = x * x' / length(x);
    peakValue = max(x .* conj(x));
    paprSymboldB(ii) = 10*log10(peakValue / meanSquareValue);
    
end

