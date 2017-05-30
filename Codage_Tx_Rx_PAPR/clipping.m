%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------ Clipping method -------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Sg_OFDM_clp = clipping(Sg_OFDM2, A_clip)

Sg_OFDM_clp = zeros(1,length(Sg_OFDM2));
	
for i = 1:length(Sg_OFDM2);
    
    if (abs(Sg_OFDM2(i)) <= A_clip);
        
        Sg_OFDM_clp(i) = Sg_OFDM2(i);
    
    else
        Sg_OFDM_clp(i) = A_clip * exp(1j*phase(Sg_OFDM2(i)));
    
    end
    
end
