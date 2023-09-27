%% Select parameters.
model_filename = 'recon2.2';  % For example: 'recon2.2', 'recon2.2withoutOR', 'mouse_recon'
expression_filename = 'NCI60_data';  % For example: 'TCGA_BRCA_data', 'NCI60_data'
data_label = 'NCI60withZielinskiBounds'; % For example: 'NCI60', 'NCI60withZielinskiBounds'
additional_constraints_filename = ''; % For example: 'Jain_data', 'NCI60_flux_bounds'
gene_set_expression_type = 'multiple'; % For example: 'multiple', 'differential'
flux_analysis_method = 'regularizedFBA'; % For example: 'FBA', 'regularizedFBA', 'pFBAnorm2'
gene_expression_map = 'exp'; % For example: 'log', 'exp'
correlation_type = 'Pearson';
experimental_reference = ''; % For example: 'proliferation', 'survival'
g = 1;
delta = 1;






%% Initialize data
load([model_filename,'.mat']);
load([expression_filename,'.mat']);

if strcmp(gene_set_expression_type, 'multiple') || strcmp(gene_set_expression_type, 'differential')
    load('reaction_expression.mat');
    load('pos_genes_in_react_expr.mat');
    load('ixs_geni_sorted_by_length.mat');
end

% if there are additional constraints, match conditions of expression data
% and of other data
if ~strcmp(additional_constraints_filename, '')
    load([additional_constraints_filename,'.mat']);
    
    common_cell_lines = intersect(data.cell_lines, bounds.cell_lines);
    [~, b] = ismember(common_cell_lines, data.cell_lines);
    expr_profile = data.feature_matrix(:, b);
    expr_genes = data.genes;
    if strcmp(experimental_reference, 'proliferation')
        proliferations = data.proliferations(b);
        tissues = data.tissues(b);
    elseif strcmp(experimental_reference, 'survival')
        survivals = data.days_to_death(data.tumour_samples_idx)'; %!
    end
    
    [~, b] = ismember(common_cell_lines, bounds.cell_lines);
    upper_bounds = bounds.upper_bounds(:, b);
    upper_bounds(upper_bounds < 0) = 0;
    lower_bounds = bounds.lower_bounds(:, b);
    lower_bounds(lower_bounds > 0) = 0;
    
else
    expr_profile = data.feature_matrix;
    expr_genes = data.genes;
    if strcmp(experimental_reference, 'proliferation')
        proliferations = data.proliferations;
        tissues = data.tissues;
    elseif strcmp(experimental_reference, 'survival')
        survivals = data.days_to_death(data.tumour_samples_idx)';
    end
    
end

[~, num_conditions] = size(expr_profile);

if strcmp(flux_analysis_method, 'FBA')
    changeCobraSolver('gurobi','LP');
elseif strcmp(flux_analysis_method, 'regularizedFBA')
    changeCobraSolver('gurobi','QP');
end

warning off





%% Get gene set expression values for all conditions
load(['C:\Users\U0033207\Desktop\MATLAB\Metabolism\Scripts\My analysis\Sensitivity analysis\num_reaction_expression_', data_label]);





%% Calculate fluxes
% set objective
objective = strcmp(model.rxns, 'biomass_reaction');
model.c(objective) = 1;
% set fluxes array
if strcmp(flux_analysis_method, 'FBA') || strcmp(flux_analysis_method, 'regularizedFBA')
    numRxns = length(model.rxns);
    fba_fluxes = zeros(num_conditions, numRxns);
elseif strcmp(flux_analysis_method, 'pFBAnorm2')
    modelIrrev = convertToIrreversible(model);
    numRxns = length(modelIrrev.rxns);
    fba_fluxes = zeros(num_conditions, numRxns);
end
feasible_conditions = true(num_conditions, 1);

for j = 1:num_conditions
    disp(j)
    target_fold_change = (num_reaction_expression(j, : ) ./ geneSet_median)';
    nan_idx = isnan(target_fold_change);
    target_fold_change(nan_idx) = 1;
    
    gamma = g * ones(length(target_fold_change), 1);
    cancer_model = model;
    %cancer_model.lb(rxns_controlled_by_genes) = cancer_model.lb(rxns_controlled_by_genes) * delta;
    %cancer_model.ub(rxns_controlled_by_genes) = cancer_model.ub(rxns_controlled_by_genes) * delta;
    if strcmp(gene_expression_map, 'log')
    for k = 1:length(target_fold_change)   %loop over the array of the geneset expressions
        if target_fold_change(k) >= 1
            cancer_model.lb(k) = cancer_model.lb(k) * (1+gamma(k)*log(target_fold_change(k)));
            cancer_model.ub(k) = cancer_model.ub(k) * (1+gamma(k)*log(target_fold_change(k)));
        else
            cancer_model.lb(k) = cancer_model.lb(k) / (1+gamma(k)*abs(log(target_fold_change(k))));
            cancer_model.ub(k) = cancer_model.ub(k) / (1+gamma(k)*abs(log(target_fold_change(k))));
        end
    end
    elseif strcmp(gene_expression_map, 'exp')
        cancer_model.lb = cancer_model.lb.*(target_fold_change.^gamma);
        cancer_model.ub = cancer_model.ub.*(target_fold_change.^gamma);
    end

    if ~strcmp(additional_constraints_filename, '')
%                 a = lower_bounds(:,j) < -0.0001;
%                 cancer_model = changeRxnBounds(cancer_model, bounds.reactions, lower_bounds(:, j), 'l');
        cancer_model = changeRxnBounds(cancer_model, bounds.reactions, upper_bounds(:,j), 'u');
    end

    % calculate fluxes 
    if strcmp(flux_analysis_method, 'FBA')
        FBAsolution = optimizeCbModel(cancer_model);
        fba_fluxes(j, :) = FBAsolution.x';
    elseif strcmp(flux_analysis_method, 'regularizedFBA')
        FBAsolution = optimizeCbModel(cancer_model, 'max', 1e-6);
        fba_fluxes(j, :) = FBAsolution.x';
    elseif strcmp(flux_analysis_method, 'pFBAnorm2')
        FBAsolution = pFBA_norm2(cancer_model, 0);
        if FBAsolution.stat == 1
            fba_fluxes(j, :) = FBAsolution.full';
        else
            feasible_conditions(j) = false;
            disp('infeasible')
        end
    end
end





%% Calculate correlation between proliferation and reactions
correlations = zeros(numRxns, 1);
pvalues = zeros(numRxns, 1);
if strcmp(experimental_reference, 'proliferation')
    for i = 1:numRxns
        [correlations(i), pvalues(i)] = corr(fba_fluxes(feasible_conditions,i), proliferations(feasible_conditions), 'type', correlation_type, 'rows', 'complete'); 
    end
elseif strcmp(experimental_reference, 'survival')
    survival = data.days_to_death(data.tumour_samples_idx)';
    for i = 1:numRxns
        [correlations(i), pvalues(i)] = corr(fba_fluxes(feasible_conditions,i), survival(feasible_conditions), 'type', correlation_type, 'rows', 'complete'); 
    end
end






%% Check highest correlations
a = pvalues <= 0.01;
high_correlations = correlations(a);
high_correlations_pvalues = pvalues(a);
if strcmp(flux_analysis_method, 'FBA') || strcmp(flux_analysis_method, 'regularizedFBA')
    correlated_reactions = model.rxns(a);
    correlated_reactions_names = model.rxnNames(a);
    correlated_subsystems = model.subSystems(a);
elseif strcmp(flux_analysis_method, 'pFBAnorm2')
    correlated_reactions = modelIrrev.rxns(a);
    correlated_reactions_names = modelIrrev.rxnNames(a);
    correlated_subsystems = modelIrrev.subSystems(a);
end
unique_subsystems = unique(correlated_subsystems);

% sort all results
[sorted_subsystems,idx] = sort(correlated_subsystems);
sorted_reactions = correlated_reactions(idx);
sorted_reactions_names = correlated_reactions_names(idx);
sorted_correlations = high_correlations(idx);
sorted_pvalues = high_correlations_pvalues(idx);

% big_table = [sorted_subsystems sorted_reactions sorted_reactions_names cell(sorted_correlations) cell(sorted_pvalues)];

fileID = fopen(['C:\Users\U0033207\Desktop\MATLAB\Metabolism\Scripts\My analysis\Sensitivity analysis\Results\', data_label, '_big_table'], 'w');
for i = 1:length(sorted_subsystems)
    fprintf(fileID,'%s\t%s\t%s\t%.2f\t%.2E\n',sorted_subsystems{i},sorted_reactions{i},sorted_reactions_names{i},round(sorted_correlations(i),2),sorted_pvalues(i));
end
fclose(fileID);

% pathway centred results
y = zeros(length(unique_subsystems),1);
z = zeros(length(unique_subsystems),1);
for i = 1:length(unique_subsystems)
    y(i) = nnz(strcmp(correlated_subsystems, unique_subsystems{i}));
    z(i) = nnz(strcmp(model.subSystems, unique_subsystems{i}));
end

% small_table = [unique_subsystems num2cell(y) num2cell(round((y.*100)./z, 2))];

fileID = fopen(['C:\Users\U0033207\Desktop\MATLAB\Metabolism\Scripts\My analysis\Sensitivity analysis\Results\', data_label, '_small_table'], 'w');
for i = 1:length(unique_subsystems)
    fprintf(fileID,'%s\t%i\t%.2f\t\n',unique_subsystems{i},y(i),round((y(i)*100)/z(i), 2));
end
fclose(fileID);

for i = 1:length(unique_subsystems)
    disp([num2str(y(i)), '                     ', unique_subsystems{i}])
end

save(['C:\Users\U0033207\Desktop\MATLAB\Metabolism\Scripts\My analysis\Sensitivity analysis\Results\', data_label, '_correlation_data4tables'], 'sorted_subsystems','sorted_reactions','sorted_reactions_names','sorted_correlations','sorted_pvalues');
