changeCobraSolver('gurobi','all');
changeCobraSolverParams('LP', 'feasTol', 1e-6);
changeCobraSolverParams('QP', 'feasTol', 1e-6);

%% load different type of models
o_auto_model = readCbModel('Data/pciCre1355/NDLadpraw_Autotrophic_Rep1.xml');
o_mixo_model = readCbModel('Data/pciCre1355/NDLadpraw_Mixotrophic_Rep1.xml');
o_hetero_model = readCbModel('Data/pciCre1355/NDLadpraw_Heterotrophic_Rep1.xml');

auto_model = o_auto_model;
mixo_model = o_mixo_model;
hetero_model = o_hetero_model;

% find exchange and uptake reactions
ex_rxns = auto_model.rxns(findExcRxns(auto_model));

% target only uptake reactions
upt_rxns = {};

for i = 1:numel(ex_rxns)

    col = find(ismember(auto_model.rxns,ex_rxns{i}));

    % also remove the demand reactions from the list
    if sum(auto_model.S(:,col)) == 1 & ~contains(ex_rxns{i}, 'DM') & contains(ex_rxns{i}, 'EX')

        upt_rxns = [upt_rxns; ex_rxns(i)];

    elseif sum(auto_model.S(:,col)) > 1

        fprintf('Need to check the stoichiometric matrix for col: %d\n',col)
        return
    end
end

% find the index for uptake reactions
auto_upt_id = find(ismember(auto_model.rxns, upt_rxns));
hetero_upt_id = find(ismember(hetero_model.rxns, upt_rxns));
mixo_upt_id = find(ismember(mixo_model.rxns, upt_rxns));

% remove the upper bound constraint of uptake reactions
for i = 1:numel(upt_rxns)
    if auto_model.ub(auto_upt_id(i)) ~= 0
        auto_model.ub(auto_upt_id(i)) = 1000;
    end

    if hetero_model.ub(hetero_upt_id(i)) ~= 0
        hetero_model.ub(hetero_upt_id(i)) = 1000;
    end

    if mixo_model.ub(mixo_upt_id(i)) ~= 0
        mixo_model.ub(mixo_upt_id(i)) = 1000;
    end
end

% find the metabolic reactions
filterRxns = @(model) model.rxns(~contains(model.rxns, {'EX_','draw','prot_pool'}));

auto_rxns   = filterRxns(auto_model);
mixo_rxns   = filterRxns(mixo_model);
hetero_rxns = filterRxns(hetero_model);

rxn_list = intersect(intersect(auto_rxns, mixo_rxns), hetero_rxns);

% remove the block reaction in rxn_list
auto_FVA = readtable('Results/FVA/auto_FVA_10p.csv');
mixo_FVA = readtable('Results/FVA/mixo_FVA_10p.csv');
hetero_FVA = readtable('Results/FVA/hetero_FVA_10p.csv');

% identify blocked reactions in the models
getBlocked = @(FVA) FVA.RxnID(FVA.minFlux == 0 & FVA.maxFlux == 0);
blocked_all = unique([getBlocked(auto_FVA); getBlocked(mixo_FVA); getBlocked(hetero_FVA)]);

rxn_list = setdiff(rxn_list, blocked_all);

% create reaction index cell array
rxn_cell = arrayfun(@(x) sprintf('Rxn%d', x), (1:numel(rxn_list))', 'UniformOutput', false);

% filter out demand and sudo reactions (containing 'No' or 'DM_')
valid_idx = ~contains(rxn_list, {'No','DM_'});
rxn_list = rxn_list(valid_idx);
rxn_cell = rxn_cell(valid_idx);

% find reversible reactions (forward and backward)
is_rev = contains(rxn_list, '_REV');
rerxn  = rxn_list(is_rev);
rerxn_cell = rxn_cell(is_rev);

forward_rxns = strrep(rerxn, '_REV', '');
forward_cells = rxn_cell(ismember(rxn_list, forward_rxns));

% combine the lists
cell_bine = [forward_cells; rerxn_cell];
rxn_bine = [forward_rxns; rerxn];

%% modeling the treatment
mutant_table = readtable('Data/Mutant_phenotypes_table_filtered_final.csv');
GO_table = readtable('Data/GO_table_filtered.txt');
% Remove the rows that contain 'antiporter'
rowsToRemove = contains(GO_table.GO_Name, 'antiporter');
GO_table(rowsToRemove, :) = [];

enzymeTable = [];

% Loop over each row in the table
for i = 1:height(mutant_table)
    % Split the 'UniProtID' by 'or'
    IDs = split(mutant_table.UniProtID{i}, ' or ');

    % Duplicate the row for each split ID
    for j = 1:numel(IDs)
        enzymeTable = [enzymeTable; IDs(j)];  % Append the new row
    end
end

go_enzyme = intersect(GO_table.UniProtID, enzymeTable);

f_mutant_table = mutant_table(contains(mutant_table.UniProtID, go_enzyme),:);

alpha = 0.1;

% change the nutrient uptake
% model1 -> auto_model
% model2 -> mixo_model
% model3 -> hetero_model
[model1, model2, model3] = changeuptake(auto_model, mixo_model, hetero_model);

% create CO2 model
% model4 -> auto_CO2
% model5 -> mixo_CO2
[model4, model1] = create_CO2_model(model1, alpha);
[model5, model2] = create_CO2_model(model2, alpha);

% create hypo model
% hypo10 -> model6
% hypo25 -> model7
% hypo75 -> model8
[model6, model7, model8, model2] = create_hypo_model(model2, alpha);

% put all models together
models = {model1; model2; model3; model4; model5; model6; model7; model8};

%% Reversible reactions
ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);
parfevalOnAll(@maxNumCompThreads, 0, 2); 

type = {'max'};
% File to store results
outputFiles = 'Results/screens/Max_flux_screen_8_Re.csv';
rxn_headers = forward_cells;

% Create headers
header = ['EnzymeID', reshape(rxn_headers, 1, [])];

% Write headers to all files
fid = fopen(outputFiles, 'w');
fprintf(fid, '%s\n', strjoin(header, ','));
fclose(fid);

for i = 1:height(f_mutant_table)

    enzyme = f_mutant_table.UniProtID{i};

    ratios = [f_mutant_table.auto_mixo(i), f_mutant_table.auto_hetero(i), f_mutant_table.mixo_hetero(i), ...
        f_mutant_table.auto_auto_CO2(i), f_mutant_table.mixo_mixo_CO2(i), f_mutant_table.mixo_hypo10_mixo(i), ...
        f_mutant_table.mixo_hypo25_mixo(i), f_mutant_table.mixo_hypo75_mixo(i)];

    % reverse log2 transform
    ratios = 2.^ratios;

    % reformat the ratios 
    ratios(6) = 1/ratios(6);
    ratios(7) = 1/ratios(7);    
    ratios(8) = 1/ratios(8);

    indices = ~isnan(ratios);

    % implement the GEM-PERSPECT
    if sum(indices) == 8
        flux_screens_8 = GEM_PERSPECT_reversible_rxns(models, ratios, rerxn, alpha, type);
        output = {enzyme, flux_screens_8{:}};
        writeToFile(outputFiles, output);
    end
end

%% Irreversible reactions
ncpu = 20;
delete(gcp('nocreate'));
parpool(ncpu);
parfevalOnAll(@maxNumCompThreads, 0, 2);

% remove the reversible from the list
rxn_list_ir = rxn_list(~ismember(rxn_list, rxn_bine));
rxn_cell_ir = rxn_cell(~ismember(rxn_cell, cell_bine));

type = {'max'};
% File to store results
outputFiles = 'Results/screens/Max_flux_screen_8.csv';
rxn_headers = rxn_cell_ir;

% Create headers
header = ['EnzymeID', reshape(rxn_headers, 1, [])];

% Write headers to all files
fid = fopen(outputFiles, 'w');
fprintf(fid, '%s\n', strjoin(header, ','));
fclose(fid);

for i = 1:height(f_mutant_table)

    enzyme = f_mutant_table.UniProtID{i};

    ratios = [f_mutant_table.auto_mixo(i), f_mutant_table.auto_hetero(i), f_mutant_table.mixo_hetero(i), ...
        f_mutant_table.auto_auto_CO2(i), f_mutant_table.mixo_mixo_CO2(i), f_mutant_table.mixo_hypo10_mixo(i), ...
        f_mutant_table.mixo_hypo25_mixo(i), f_mutant_table.mixo_hypo75_mixo(i)];

    % reverse log2 transform
    ratios = 2.^ratios;

    % reformat the ratios 
    ratios(6) = 1/ratios(6);
    ratios(7) = 1/ratios(7);    
    ratios(8) = 1/ratios(8);

    indices = ~isnan(ratios);

    % implement the GEM-PERSPECT
    if sum(indices) == 8
        flux_screens_8 = GEM_PERSPECT(models, ratios, rxn_list_ir, alpha, type);
        output = {enzyme, flux_screens_8{:}};
        writeToFile(outputFiles, output);
    end
end
