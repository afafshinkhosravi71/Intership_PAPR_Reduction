%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------- Calcul PAPR ------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PaprSymb_dB = calculPAPR(Sg_final2, NSymb)
 
PaprSymb_dB = zeros(1, NSymb);
% Parcourir chaque symbole 
for i = 1:1:NSymb;
    
    temp = transpose(Sg_final2(:,i));
        
    % Peak to average power ratio PAPR
    MeanSquareValue = temp * temp' / length(temp);
    PeakValue = max(temp .* conj(temp));
    PaprSymb_dB(i) = 10 * log10(PeakValue / MeanSquareValue) ;    
    
end