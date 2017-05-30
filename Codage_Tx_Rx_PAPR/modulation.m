%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------------- Modulation du flux binaire ---------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sg_Mod = modulation (NbSym_Mod, Msg_bin)

Sg_Mod = modulate(modem.qammod('M',NbSym_Mod,'InputType','Bit'),Msg_bin); 

% Constelleation
scatterplot(Sg_Mod);
grid on;
title('Constellation apres modulation');

