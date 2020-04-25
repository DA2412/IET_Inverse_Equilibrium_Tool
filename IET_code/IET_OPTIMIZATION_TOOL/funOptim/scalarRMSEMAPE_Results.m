function [MAPE,RMSE] = scalarRMSEMAPE_Results(W,G_boundary,PSI_coils_boundary,Psi_axis,Psi_boundary)

residuals = G_boundary*W' - PSI_coils_boundary;

DeltaPsi = abs(Psi_axis - Psi_boundary);
Norm_factor = DeltaPsi;

MAPE = max(abs(residuals))./Norm_factor*100;
RMSE = sqrt(mean((residuals).^2))./Norm_factor*100;
% w1 = lambda;
% w2 = 1-lambda;

% fobj = w1 * RMSE + w2 * MAPE;
% W_obj = W;
end

