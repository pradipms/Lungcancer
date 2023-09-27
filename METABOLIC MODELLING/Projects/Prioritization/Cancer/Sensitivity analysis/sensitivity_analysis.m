%% Select parameters.
model_filename = 'recon2.2';  % For example: 'recon2.2', 'recon2.2withoutOR', 'mouse_recon'
expression_filename = 'NCI60_data';  % For example: 'GTEx_data_53_tissues', 'NCI60_data'
data_label = 'NCI60withZielinskiBounds'; % For example: 'NCI60', 'NCI60withZielinskiBounds'
additional_constraints_filename = ''; % For example: 'Jain_data', 'NCI60_flux_bounds'
gene_set_expression_type = 'multiple'; % For example: 'multiple', 'differential'
flux_analysis_method = 'regularizedFBA'; % For example: 'FBA', 'regularizedFBA'
gene_expression_map = 'exp'; % For example: 'log', 'exp'
correlation_type = 'Pearson';
experimental_reference = ''; % For example: 'proliferation', 'survival'
generate_reaction_expression = true;
save_results = false;
gamma_values = [1];
% gamma_values = [1 3 10 30 100 300 1000 3000 10000 30000];
% gamma_values = [1 3 10 30 50 70 80 90 100 110 111 112 113 114 115 116 117 118 119 120 130 140 150 200 300 1000 3000 10000 30000 100000];
delta_values = [1]; % native recon bound rescaling to make uptake and gene expression constraints comparable
% delta_values = [0.000025:0.000001:0.000035];
% delta_values = [0.000001 0.000003 0.00001 0.00003 0.0001 0.0003 0.001 0.003];


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
    changeCobraSolver('ibm_cplex','LP');
elseif strcmp(flux_analysis_method, 'regularizedFBA')
    changeCobraSolver('ibm_cplex','QP');
end

warning off






%% Get gene set expression values for all conditions
if generate_reaction_expression == true
    disp('Generating gene set expressions...')
    
    to_replace = 'min.\d*+\.+\d*,\d*+\.+\d*.|max.\d*+\.+\d*,\d*+\.+\d*.|min..\d*+\.+\d*.,\d*+\.+\d*.|max..\d*+\.+\d*.,\d*+\.+\d*.|min..\d*+\.+\d*.,.\d*+\.+\d*..|max..\d*+\.+\d*.,.\d*+\.+\d*..|min.\d*+\.+\d*,.\d*+\.+\d*..|max.\d*+\.+\d*,.\d*+\.+\d*..';  %searches for all the strings of kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM) or  min((NUM.NUM),NUM.NUM) or max((NUM.NUM),NUM.NUM) or  min((NUM.NUM),(NUM.NUM)) or max(NUM.NUM,(NUM.NUM)) or  min(NUM.NUM,(NUM.NUM)) or max((NUM.NUM),(NUM.NUM))
    to_replace_nan1 = '|min.NaN,\d*+\.+\d*.|min.\d*+\.+\d*,NaN.|max.NaN,\d*+\.+\d*.|max.\d*+\.+\d*,NaN.|min..NaN.,\d*+\.+\d*.|min..\d*+\.+\d*.,NaN.|max..NaN.,\d*+\.+\d*.|max..\d*+\.+\d*.,NaN.|min..NaN.,.\d*+\.+\d*..|min..\d*+\.+\d*.,.NaN..|max..NaN.,.\d*+\.+\d*..|max..\d*+\.+\d*.,.NaN..|min.NaN,.\d*+\.+\d*..|min.\d*+\.+\d*,.NaN..|max.NaN,.\d*+\.+\d*..|max.\d*+\.+\d*,.NaN..';  %searches for all the strings of kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM) or  min((NUM.NUM),NUM.NUM) or max((NUM.NUM),NUM.NUM) or  min((NUM.NUM),(NUM.NUM)) or max(NUM.NUM,(NUM.NUM)) or  min(NUM.NUM,(NUM.NUM)) or max((NUM.NUM),(NUM.NUM))
    to_replace_nan2 = '|min.NaN,NaN.|max.NaN,NaN.|min..NaN.,NaN.|max..NaN.,NaN.|min..NaN.,.NaN..|max..NaN.,.NaN..|min.NaN,.NaN..|max.NaN,.NaN..'; %searches for all the strings of kind min(NUM.NUM,NUM.NUM) or max(NUM.NUM,NUM.NUM) or  min((NUM.NUM),NUM.NUM) or max((NUM.NUM),NUM.NUM) or  min((NUM.NUM),(NUM.NUM)) or max(NUM.NUM,(NUM.NUM)) or  min(NUM.NUM,(NUM.NUM)) or max((NUM.NUM),(NUM.NUM))
    to_replace = [to_replace to_replace_nan1 to_replace_nan2];
    
    pos_genes_in_dataset = zeros(numel(model.genes), 1);
    x = nan(num_conditions, numel(model.genes));
    for k = 1:numel(model.genes)
        position = find(strcmp(model.genes{k}, expr_genes));
        if ~isempty(position)
            pos_genes_in_dataset(k) = position;
            x(:, k) = expr_profile(pos_genes_in_dataset(k), :);
        end
    end

    num_reaction_expression = zeros(num_conditions, length(reaction_expression));

    for t = 1:num_conditions
        % Calculate expression bound for each reaction
        %indices are sorted by length of their string, so the longest get replaced first. This avoids, for example, that if we have two genes 123.1 and 23.1, the substring '23.1' is replaced in both. If we first replace the longest and finally the shortest, this problem is avoided
        reaction_expr = reaction_expression;
        for i = ixs_geni_sorted_by_length %loop over the array of the non-1 gene expressions, in order to replace the names of genes in geni_reazioni.mat with their values. All the gene set expressions with only 1s as gene values , will be left empty and at the end of this loop everything empty will be substituted with 1 anyway. This avoids looping over all the genes yt, which is very expensive
            posizioni_gene = pos_genes_in_react_expr{i};
            for j = 1:length(posizioni_gene)      %for each string of reaction_expression, we replace the substring 'bXXXX' with the number representing its gene expression
                reaction_expr{posizioni_gene(j)} = strrep(reaction_expr{posizioni_gene(j)}, model.genes{i}, num2str(x(t, i),'%.15f'));  %Matlab strangely truncates decimal digits when using num2str. Addimg %.12f at least ensures that 12 decimal digits are included in the number converted into string
            end
        end
        reaction_expr( cellfun(@isempty, reaction_expr) ) = {'1.0'};  %replaces all the empty cells of gene expression (e.g. exchange reactions) with 1, i.e. gene expressed nomally

        for i = 1:length(num_reaction_expression)
            str = reaction_expr{i};
            num_parenthesis = numel(strfind(str, ')'));
            while (num_parenthesis > 32) %if there are more than 32 parentheses, matlab is unable to run EVAL. So we need to reduce these parentheses manually by starting to eval smaller pieces of the string
                substrings_to_replace = regexp(str, to_replace, 'match');
                if isempty(substrings_to_replace)
                    num_parenthesis = 0; %if num_parenthesis > 32 and there is nothing that can be replaced with regexp, we force this, in order to avoid an endless loop. Later, eval will catch an exception as it cannot evaluate when num_parenthesis>32
                else
                    for j = 1:numel(substrings_to_replace)
                        ss_rep = substrings_to_replace{j};
                        str = strrep(str, ss_rep, num2str(eval(ss_rep), '%.15f'));
                    end
                    num_parenthesis = numel(strfind(str, ')'));
               end
            end

            str = regexprep(str, '/', '');
            num_reaction_expression(t, i) = eval(str);   %evaluates the cells like they are numerical expressions (so as to compute min and max of gene expressions)
        end

    end
    
    
    if strcmp(gene_set_expression_type, 'multiple')
        % Calculate mean of gene set expression
        geneSet_mean = mean(num_reaction_expression, 'omitnan');
        geneSet_median = median(num_reaction_expression, 'omitnan');
    
    elseif strcmp(gene_set_expression_type, 'differential')
        % Calculate mean of normal samples expression
        geneSet_mean = mean(num_reaction_expression(data.normal_samples_idx, :), 'omitnan');
        geneSet_median = median(num_reaction_expression(data.normal_samples_idx, :), 'omitnan');

    end
    
    
    save(['C:\Users\U0033207PC\Desktop\MATCANCER\Metabolism\Scripts\My analysis\Sensitivity analysis\num_reaction_expression_', data_label], ...
        'num_reaction_expression', 'geneSet_mean', 'geneSet_median');
    disp('Done')
    
else
    load(['C:\Users\U0033207PC\Desktop\MATCANCER\Projects\Prioritization\Cancer\Sensitivity analysis\num_reaction_expression_', data_label]);

end






%% Perform sensitivity analysis
if strcmp(gene_set_expression_type, 'differential')
    num_reaction_expression = num_reaction_expression(data.tumour_samples_idx, :);
    num_conditions = sum(data.tumour_samples_idx);
end

% set growth reactions
model = changeObjective(model, 'biomass_reaction');
objective_idx = strcmp(model.rxns, 'biomass_reaction');
num_gammas = length(gamma_values);
num_deltas = length(delta_values);
predicted_growths = zeros(num_conditions, num_gammas, num_deltas);
rxns_controlled_by_genes = ~cellfun(@isempty, model.grRules);

% loop over gamma values
for i = 1:num_gammas
    disp(['gamma = ', num2str(gamma_values(i))])
    
    for l = 1:num_deltas
        
        % loop over conditions
        for j = 1:num_conditions
            % generate condition-specific model
            target_fold_change = (num_reaction_expression(j, :) ./ geneSet_median)';
            nan_idx = isnan(target_fold_change);
            target_fold_change(nan_idx) = 1;

            gamma = gamma_values(i) * ones(length(target_fold_change), 1);
            cancer_model = model;
            cancer_model.lb(rxns_controlled_by_genes) = cancer_model.lb(rxns_controlled_by_genes) * delta_values(l);
            cancer_model.ub(rxns_controlled_by_genes) = cancer_model.ub(rxns_controlled_by_genes) * delta_values(l);
            
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
                cancer_model = changeRxnBounds(cancer_model, bounds.reactions, lower_bounds(:, j), 'l');
%                 cancer_model = changeRxnBounds(cancer_model, bounds.reactions, upper_bounds(:,j), 'u');
            end

            % calculate fluxes
            if strcmp(flux_analysis_method, 'FBA')
                FBAsolution = optimizeCbModel(cancer_model);
            elseif strcmp(flux_analysis_method, 'regularizedFBA')
                FBAsolution = optimizeCbModel(cancer_model, 'max', 1e-6);
            end
    %         C_model.S = [C_model.S; C_model.c'];
    %         C_model.b = [C_model.b; FBAsolution.f];
    %         C_model.c = zeros(length(C_model.c), 1);
    %         second_objective = find(~cellfun(@isempty, strfind(model.rxns, 'ATPS4m')));
    %         model.c(second_objective) = 1;
    %         FBAsolution = optimizeCbModel(C_model, 'max');

            if FBAsolution.stat == 1
                predicted_growths(j,i,l) = sum(FBAsolution.v(objective_idx));
            end
        end
        
        if sum(predicted_growths(:,i,l)) == 0
            disp('   no growth delta = ', num2str(delta_values(l)))
        end
    end
end



% evaluate predicted biomass/fluxes in relation to experimental values
correlations = zeros(num_gammas, num_deltas);
pvalues = zeros(num_gammas, num_deltas);

if strcmp(experimental_reference, 'proliferation')
    for i = 1:num_gammas
        for j = 1:num_deltas
        [correlations(i,j), pvalues(i,j)] = corr(predicted_growths(:,i,j), proliferations, 'type', correlation_type, 'tail', 'right'); 
        end
    end
    
    if save_results == true
        save(['C:\Users\U0033207PC\Desktop\MATCANCER\Metabolism\Scripts\My analysis\Sensitivity analysis\Results\sensitivity_', ...
            data_label, '_', flux_analysis_method], ...
            'predicted_growths', 'proliferations', 'correlations', 'pvalues', 'gamma_values', 'delta_values', 'tissues');
    end
    
    
elseif strcmp(experimental_reference, 'survival')
    for i = 1:num_gammas
        for j = 1:num_deltas
        [correlations(i,j), pvalues(i,j)] = corr(predicted_growths(:,i,j), survivals, 'type', correlation_type, 'tail', 'left', 'rows', 'complete'); 
        end
    end
    
    if save_results == true
        save(['C:\Users\U0033207PC\Desktop\MATCANCER\Metabolism\Scripts\My analysis\Sensitivity analysis\Results\sensitivity_', ...
            data_label, '_', flux_analysis_method], ...
            'predicted_growths', 'survivals', 'correlations', 'pvalues', 'gamma_values', 'delta_values');
    end
end






%% Plot
figure(1)
h = plot(gamma_values, correlations, 'Marker', 'o');
set(gca,'xscale','log')
xlabel('\gamma', 'FontSize', 25)
ylabel('$r$', 'FontSize', 25, 'Interpreter', 'latex')
set(gca,'LineWidth',1.5)
set(h(:), 'LineWidth', 1.5)
legend(num2str(delta_values'))

figure(2)
if exist('proliferations', 'var')
    [max_gamma_idx, max_delta_idx] = ind2sub([num_gammas num_deltas], find(correlations == max(correlations(:))));
    colormap jet
    h = gscatter(proliferations, predicted_growths(:, max_gamma_idx, max_delta_idx), tissues, [], 'os+x*^v<>', 10);
    xlabel('measured proliferation [h^{-1}]')
    ylabel('predicted biomass yield [h^{-1}]')
    set(gca,'FontSize',20,'LineWidth',1.5)
    set(h(:), 'LineWidth', 1.5)
    legend('Location','northeastoutside')
    
elseif exist('survivals', 'var')
    [min_gamma_idx, min_delta_idx] = ind2sub([num_gammas num_deltas], find(correlations == min(correlations(:))));
    colormap jet
    h = gscatter(survivals, predicted_growths(:, min_gamma_idx, min_delta_idx), [], [], 'o', 15);
    xlabel('survival [d]')
    ylabel('predicted biomass yield [h^{-1}]')
    set(gca,'FontSize',20,'LineWidth',1.5)
    set(h(:), 'LineWidth', 1.5)
    
else
    [max_gamma_idx, max_delta_idx] = ind2sub([num_gammas num_deltas], find(correlations == max(correlations(:))));
    scatter(proliferations, predicted_growths(:, max_gamma_idx, max_delta_idx))
    xlabel('measured proliferation [h^{-1}]')
    ylabel('predicted biomass yield [h^{-1}]')
    set(gca,'FontSize',15)
end
