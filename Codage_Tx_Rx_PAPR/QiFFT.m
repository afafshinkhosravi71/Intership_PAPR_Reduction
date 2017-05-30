%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------- IFFT du signal --------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sg_OFDM2 = QiFFT(Sg_OFDM, NSymb, NFFTSize)

% Application de la iFFT
% Construction de la matrice Q_ifft de Fourrier inverse (64 , 64)
Q_iFFT = zeros(NFFTSize, NFFTSize);
for i = 0:NFFTSize-1
    for k = 0:NFFTSize-1
        Q_iFFT(i+1, k+1) = exp(2*pi*1j*k*i / NFFTSize);
    end
end

Sg_OFDM2 = [];

% Parcourir chaque symbole 
for i = 1:1:NSymb;
    
    % Recupere les trames OFDM
    frame = Sg_OFDM(:,i);   
    
    % Application de la iFFT de nFFFTSize 
    Sg_IFFT = transpose(Q_iFFT * frame);   
    
    %% Domaine temporel 
    % Concatenation des symboles
    Sg_OFDM2 = [Sg_OFDM2 Sg_IFFT];   
    
end