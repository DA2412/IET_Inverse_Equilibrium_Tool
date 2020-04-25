function [fobj,W] = onlymu1(W_norm,v_lb,v_ub,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary)
% load('./output_data/OPTIM_DATA_INPUT.mat', 'G_boundary','PSI_coils_boundary',...
%     'psiBoundary','psi_ak');

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

residuals = bsxfun(@minus,G_boundary*W,PSI_coils_boundary);

DeltaPsi = abs(Psi_axis - Psi_boundary);
Norm_factor = DeltaPsi;


f1 = max(abs(residuals),[],1)./Norm_factor*100;

if l_flip
    fobj =  f1';
else
    fobj = f1;
end

end

