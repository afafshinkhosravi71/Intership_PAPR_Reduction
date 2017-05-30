function [papr_dB, meanSquare, peakVal] = Calcul_papr(input_OFDM, nSymbol)

nBitperSymbol = 80;  % nFFTSize(64) + n_IG(16) = 80

input_OFDM = reshape(input_OFDM, nBitperSymbol, nSymbol); 

for ii = 1:nSymbol
    
    x = transpose(input_OFDM(:,ii));
        
    % Calcul du peak to average power ratio PAPR
    meanSquare(ii) = x * x' / length(x);
    peakVal(ii) = max(x .* conj(x));
    papr_dB(ii) = 10*log10(peakVal(ii) / meanSquare(ii));
    
end

