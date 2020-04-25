% Here comes the non-linear constraints.
function [C, Ceq]=SQPnolcon(W_norm,Utopia,L,N,Xpj,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary)


mubarrato=muNorm(W_norm,Utopia,L,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary);

     
        %C=N * (mubarrato-Xpj)';%This one is required in NNC.
        C=N' * (mubarrato-Xpj);%This one is required in NNC.

        Ceq=0;
  
end