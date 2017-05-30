%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---  Programme de calcul des fonctions de sensibilites des parametres de l'amplificateur  ---- 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SENS] = sensb(ysim, A_clip, L, nb_para, nFFTSize) 

indice_sup = find(abs(ysim) > A_clip); 
eps = A_clip - abs(ysim(indice_sup)); 
Nb_point = length(indice_sup);

F = zeros(Nb_point, nFFTSize);

for i = 1:Nb_point;    
    F(i, indice_sup(i)) = exp(1j*phase(ysim(indice_sup(i))));    
end;

SENS = (F * L);

% SENS = [];
% 
% j = 1;
% for m = 1:length(ysim)  
%     if (abs(ysim(m)) > A_clip)
%         for i = 1:nb_para
%             SENS(j,i) =  L(m,i) * exp(1j*phase(ysim(m)));   %  Calcul du Gradient
%         end
%         j = j + 1;
%     end;
% end;
