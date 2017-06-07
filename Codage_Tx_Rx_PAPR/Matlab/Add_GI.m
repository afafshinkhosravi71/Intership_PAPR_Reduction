%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------ Ajout d'intevalle de garde 16 symboles ---------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sg_final2, Nb_GI] = Add_GI(Sg_OFDM2, NFFTSize, NSymb)

% Dispose les données en matrice (NFFTSize, NSymb) (64 , 22)
Sg_OFDM2 = reshape(Sg_OFDM2, NFFTSize, NSymb); 

% Nombre de sous porteuses pour l'intervalle de garde 
% Pourcentage = 100 * Nb/GI / nFFTSize;
Nb_GI = 16;    

% ================================================================ %
% Methode 1 : Copie des 16 deniers symboles dans l'ig (time)
% ================================================================ %
Sg_final2 = [];

% Parcourir chaque symbole 
for i = 1:1:NSymb;
    
    temp = Sg_OFDM2(:,i);
    
    GI = transpose(temp(49:64));   % last samples  
    
    % Rajoute des prefixes supplementaires de 16 echantillons  
    frame_GI = [GI  transpose(temp)];  
    
    % Concatenation des symboles pour former la sortie finale
    Sg_final2 = [Sg_final2 frame_GI];
    
end