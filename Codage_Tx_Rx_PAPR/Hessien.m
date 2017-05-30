%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--- Simulation avec la methode TR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [input_iFFT_TR] = Hessien(input_iFFT_mat, A_clip, iter_max,nSymbol,Gain_OFDM)

deb = 1; fin = nSymbol;
nFFTSize = 64;            %% debit= 6Mbits/s;

% ------------------------------------------------------------------------
%    Identification des parametres C avec l'algorithme gradient conjugu� (direction 1) 
% -----------------------------------------------------------------------

%%-------- Construction de la matrice Q_ifft de Fourrier inverse
for i = 0:nFFTSize-1 ; 
    for k = 0:nFFTSize-1 ; 
        Q_iFFT(i+1, k+1) = exp(2*pi*1j*k*i/nFFTSize) ; 
    end; 
end;

%%----------- Disposition des porteuses reservees a la TR
% 12 porteuses reservees a la TR [ 6 -  1  - 5 ]
%TR_porteuse = ones(64,1);
TR_porteuse = ([ones(6,1); zeros(26,1); ones(1,1);zeros(26,1); ones(5,1)]); 

%=================================================================%
%%------- PAPR_Debit_6 Mbits/s;
%=================================================================%
input_iFFT_TR = []; 

HH = diag(TR_porteuse); H = []; 
for i=1:length(HH); 
    if (sum(HH(:,i)) == 1); 
        H = [H HH(:,i)]; 
    end
end

[nb_row nb_col] = size(H); 
nb_para = nb_col;

L = Q_iFFT * H;

%===========================
%%----- Filtre, pr�fixes et sufixes 
%==========================

for ii = deb:fin;
    
    Crit = []; Theta = [];  
    
    inputiFFT = input_iFFT_mat(:,ii) / Gain_OFDM;  
    % Transpose les sous porteuses des indices [-26 to -1] aux indices [38 to 63]
    
    x = Q_iFFT * (inputiFFT) ; 
%     disp('taille de x')
%     size(x)
    % M�thode 2 : iFFT � partir de la fonction iFFFT de Matlab 
    
    % ---------------------------------------------------------------------------------------------
    %  Recherche des valeurs initiales
    % ---------------------------------------------------------------------------------------------
    
    y_clip = x; 
    indice_sup = find(abs(y_clip) > A_clip); 
    y_clip(indice_sup) = exp(1j*phase(y_clip(indice_sup))) .* A_clip;
    
    iter  = 1;   
    Theta(:,iter) = ((L'*L)^(-1)*L'*(y_clip - x));
    %Theta(:,iter) = zeros(nb_para,1);
    
    ysim = x + Q_iFFT * (H * Theta(:,iter));
    
    % ---------------------------------------------------------------------------------------------
    %  E S T I M A T I O N    P A R A M E T R I Q U E   P N L :  I N I T I A L I S A T I O N S
    % ---------------------------------------------------------------------------------------------
          
    % ------------------------------------------------------------------------------
    %                  ALGORITHME DE MARQUARDT
    % ------------------------------------------------------------------------------

    % ------ Calcul du crit�re initial ---------
   
    indice_sup = find(abs(ysim) > A_clip);    
    eps = A_clip - abs(ysim(indice_sup)); Nb_point = length(indice_sup);
    Crit(iter) = eps' * eps ;
    
    test = Nb_point > 0;

    while (test),

     % ------------ Calcul du gradient et du Hessien � l'it�ration n�iter ---------------

       [SENS] = sensb(ysim, A_clip, L, nb_para, nFFTSize) ;    % Matrice des fonctions de sensibilit�

       Grad = -(transpose(SENS) * eps)  ;            %  Calcul du Gradient       
       Hess  = (transpose(L) * L) ;              %  Calcul du Hessien 

         % actualisation de l'�cart delta       

          delta = -inv(Hess) * Grad ;
          Theta_recherche = Theta(:,iter) + delta ;

        % ------- Calcul des nouveaux �carts -----------
        
          inputiFFT_TR = inputiFFT + (H * Theta_recherche);
          ysim = Q_iFFT * (inputiFFT_TR); %           
          indice_sup = find(abs(ysim) > A_clip); eps = A_clip - abs(ysim(indice_sup));
          
          % -------- Calcul du crit�re ----------

             Critere = eps' * eps;             
       
       test = (iter <= iter_max && ~isempty(indice_sup));
       
       iter=iter+1; 

       Theta(:,iter) = Theta_recherche;           % actualisation de Theta
       Crit(iter) = Critere ; 

    end

    % ------------------------------------------------------------------------------
    %                  FIN MARQUARDT
    % ---------------------------------------------------------------------
    
    
    % Concatenation des symboles pour former la sortie finale
    input_iFFT_TR = [input_iFFT_TR inputiFFT_TR];     
    
end

% Visualisation
figure()
subplot(211)
stem(real(input_iFFT_mat(1:64,1)),'r*')
title('Trame d''un symbole OFDM avant TR')
subplot(212)
stem(real(input_iFFT_TR(1:64)), '*')
title('Trame d''un symbole OFDM apres TR')


%plot(Crit)

%hold on; plot(abs(x),'b'); plot(abs(y_clip),'r');
%hold on; plot(abs(x),'b'); plot(abs(ysim),'r');


%input_iFFT_TR = reshape(input_iFFT_TR,nFFTSize,nSymbol); 