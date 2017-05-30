%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------- Modulation OFDM IEEE 802.11a -------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Sg_OFDM, NSymb ] = Allocation_OFDM(Sg_Mod, NFFTSize)

% Taille des données modulées
NSymbol_Mod = length(Sg_Mod);

% (1 symbol OFDM 64 = 52 Subcarriers(48  data + 4 pilots) + 12 GI a zero (intervalle de garde)

NSubcarriersData = 48;       % Nombre de sous-porteuse OFDM pour les data
IndexSubcarrierData = [-26:-22 -20:-8 -6:-1 1:6 8:20 22:26];    % Emplacements des 52 sous-porteuses data
IndexSubcarrierPilot = [-21 -7 7 21];                           % Emplacements des 4 sous-porteuses pilote

NSymb = ceil(NSymbol_Mod/NSubcarriersData);                     % nombre de symbole arrondi au superieur 

% Complete les bits manquants par zero pour atteindre NSubcarriersData * NSymb
Sg_Mod2 = [Sg_Mod ; zeros(NSubcarriersData * NSymb - NSymbol_Mod ,1)];  

% Dispose les données en matrice (NSymb , NSubcarriersData) pour la transmission OFDM
Sg_Mod2 = reshape(Sg_Mod2, NSymb, NSubcarriersData);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%----------------- Insertion des données et des pilotes ------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sg_OFDM = []; 
Pilot_Value = [];
% Parcourir chaque symbole 
for i = 1:1:NSymb;
    
    % Initialisation du symbole (1, 64)
    Symb = zeros(1, NFFTSize);  
    
    % Assignation des 48 bits des datas sur les subcarriers_data 48
    Symb(1, IndexSubcarrierData + NFFTSize/2 + 1) = Sg_Mod2(i,:);

    % Ajout des piloteSg_OFDMs sur les sous-porteuses 4
    % Generation aletoire de pilote (complex)
    Pilot = (((rand(1,4)' > 0.5) - 0.5) + 1j*((rand(1,4)' > 0.5) - 0.5))*2;   
    
    Symb(1, IndexSubcarrierPilot + NFFTSize/2 + 1) = Pilot;    
    
    % Concatenation des symboles pour former la sortie finale
    Sg_OFDM = [Sg_OFDM Symb];  
    
    % Recuperation des pilotes
    Pilot_Value = [Pilot_Value Pilot];
    
end

% Visualisation
figure()
stem(real(Sg_OFDM(1:64)), '*')
title('Trame d''un symbole OFDM')

% Dispose les données en matrice (NFFTSize, NSymb) (64 , 22)
Sg_OFDM = reshape(Sg_OFDM, NFFTSize, NSymb); 