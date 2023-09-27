% MATLAB 2017b
% COBRA toolbox version 3.1
% IBM CPLEX version 12.8

load recon2.2.mat
changeCobraSolver('ibm_cplex', 'all');
model = changeObjective(model, 'biomass_reaction');
solution = optimizeCbModel(model, 'max', 1e-6);
fluxes = solution.v;
xlswrite('recon_control_fluxes', fluxes);