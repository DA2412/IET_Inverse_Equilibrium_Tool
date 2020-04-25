function [fobj] = mu2Norm(W_norm,Utopia,L,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary)

% Questa parte l'ho messa perche' a quanto pare fmincon e GA prendono 
% gli input in modo diverso e mi serviva una funzione sola per fare dei 
% test
v_size = size(W_norm);
if v_size(1) ~= numel(v_lb)
    W_norm = W_norm';
    l_flip = true;
else
    l_flip = false;
end

W = bsxfun(@plus,v_lb',bsxfun(@times,W_norm,(v_ub-v_lb)'));

residuals = G_boundary*W - PSI_coils_boundary;

DeltaPsi = abs(Psi_axis - Psi_boundary);
Norm_factor = DeltaPsi;

mu1 = max(abs(residuals))/Norm_factor*100;
mu2 = sqrt(mean((residuals).^2))/Norm_factor*100;

f = [(mu1-Utopia(1))/L(1); (mu2 - Utopia(2))/L(2)];
fobj = f(2);


end

