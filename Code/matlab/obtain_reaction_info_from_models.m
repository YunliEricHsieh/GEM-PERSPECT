% load models
o_auto_model = readCbModel('Data/pciCre1355/NDLadpraw_Autotrophic_Rep1.xml');
o_mixo_model = readCbModel('Data/pciCre1355/NDLadpraw_Mixotrophic_Rep1.xml');
o_hetero_model = readCbModel('Data/pciCre1355/NDLadpraw_Heterotrophic_Rep1.xml');

auto_model = o_auto_model;
mixo_model = o_mixo_model;
hetero_model = o_hetero_model;

% find the metabolic reactions
filterRxns = @(model) model.rxns(~contains(model.rxns, {'EX_','draw','prot_pool'}));

auto_rxns   = filterRxns(auto_model);
mixo_rxns   = filterRxns(mixo_model);
hetero_rxns = filterRxns(hetero_model);

rxn_list = intersect(intersect(auto_rxns, mixo_rxns), hetero_rxns);

% read FVA data
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

%% find reversible reactions (forward and backward)
is_rev = contains(rxn_list, '_REV');
rerxn  = rxn_list(is_rev);
rerxn_cell = rxn_cell(is_rev);

forward_rxns = strrep(rerxn, '_REV', '');
forward_cells = rxn_cell(ismember(rxn_list, forward_rxns));

% combine the lists
cell_bine = [forward_cells; rerxn_cell];
rxn_bine = [forward_rxns; rerxn];

% save table
final_table = cell2table([cell_bine, rxn_bine], 'VariableNames', {'RxnIndex', 'RxnID'});
writetable(final_table, 'Data/Reactions/list_of_reversible_rxns.csv');

%% find irreversible reactions
% remove the reversible from the list
rxn_list_ir = rxn_list(~ismember(rxn_list, rxn_bine));
rxn_cell_ir = rxn_cell(~ismember(rxn_cell, cell_bine));

% save table
final_table = cell2table([rxn_cell_ir, rxn_list_ir], 'VariableNames', {'RxnIndex', 'RxnID'});
writetable(final_table, 'Data/Reactions/list_of_irreversible_rxns.csv');

%% find reactions without GPR rules
% find indices of reactions from rxn_list in  model
[~, a_idx] = ismember(rxn_list, auto_model.rxns);

% identify reactions without GPR rules (empty rules) in model
auto_noGPR   = cellfun(@isempty, auto_model.rules(a_idx));

% remove the backward reactions 
no_REV = ~contains(rxn_list, '_REV');

% combine both criteria (no GPR rules AND no '_REV')
selected_rxns = auto_noGPR & no_REV;

% save table
auto_table   = cell2table([rxn_cell(selected_rxns), rxn_list(selected_rxns)], 'VariableNames', {'RxnIndex', 'RxnID'});
writetable(auto_table, 'Data/Reactions/list_of_rxns_without_GPR.csv');

%% find the rxns only associated with the enzyme that existing in models and mutant table
mutant_table = readtable('Data/condition_table_filtered.csv');
GO_table = readtable('Data/GO_table_filtered.txt');

% remove the rows that contain 'antiporter'
rowsToRemove = contains(GO_table.GO_Name, 'antiporter');
GO_table(rowsToRemove, :) = [];

enzymeTable = [];

% loop over each row in the table
for i = 1:height(mutant_table)
    % split the 'UniProtID' by 'or'
    IDs = split(mutant_table.UniProtID{i}, ' or ');

    % duplicate the row for each split ID
    for j = 1:numel(IDs)
        enzymeTable = [enzymeTable; IDs(j)];  % Append the new row
    end
end

go_enzyme = intersect(GO_table.UniProtID, enzymeTable);

% filter enzymes in models
pattern = {'unknown', 'CHLRE', 'ChreC'};
model_genes = intersect(intersect(auto_model.genes(~contains(auto_model.genes, pattern)), ...
                       hetero_model.genes(~contains(hetero_model.genes, pattern))), ...
                       mixo_model.genes(~contains(mixo_model.genes, pattern)));

model_enzyme_in_table = intersect(go_enzyme, model_genes);
target_number = find(ismember(auto_model.genes, model_enzyme_in_table));

% identify reactions with only target enzymes
isValid = false(numel(auto_model.rules), 1);
for i = 1:numel(auto_model.rules)
    numbers = regexp(auto_model.rules{i}, 'x\((\d+)\)', 'tokens');
    if ~isempty(numbers)
        nums = str2double([numbers{:}]);
        isValid(i) = all(ismember(nums, target_number));
    end
end

target_rxns = intersect(auto_model.rxns(isValid), rxn_list);

% remove backward reactions
target_rxns = target_rxns(~contains(target_rxns, '_REV'));

% find the genes in GPR rules of target rxns
enzyme_list = [];

for i = 1:numel(target_rxns)
    ind = find(ismember(auto_model.rxns, target_rxns{i}));
    t_enzyme = regexp(auto_model.rules{ind}, 'x\((\d+)\)', 'tokens');
    numbers = cellfun(@(x) str2double(x{1}), t_enzyme);
    tmp_cell = [];
    for j = 1:numel(numbers)
        tmp_cell = [tmp_cell, auto_model.genes(numbers(j))];
    end
    enzyme_list{i} = unique(tmp_cell);
end

enzyme_list = enzyme_list';

for i = 1:numel(enzyme_list)
    if iscell(enzyme_list{i})
        enzyme_list{i} = strjoin(enzyme_list{i}, ';');
    end
end

% find the rxns index
rxnIndex = rxn_cell(ismember(rxn_list, target_rxns));

% Create final table and write to file
final_table = cell2table([rxnIndex, target_rxns, enzyme_list], ...
    'VariableNames', {'RxnIndex', 'RxnID', 'Enzymes'});

writetable(final_table, 'Data/Reactions/list_of_rxns_with_enzymes_in_both.csv');

%% find the rxns associated with at least one enzyme that existing in model and mutant table
% identify reactions with at least one target enzymes
isValid = false(numel(auto_model.rules), 1);
for i = 1:numel(auto_model.rules)
    matches = regexp(auto_model.rules{i}, 'x\((\d+)\)', 'tokens');
    if ~isempty(matches)
        numbers = cellfun(@(x) str2double(x{1}), matches);
        % Check if at least one number is in target list
        if any(ismember(numbers, target_number))
            isValid(i) = true;
        end
    end
end

target_rxns = intersect(auto_model.rxns(isValid), rxn_list);

% remove backward reactions
target_rxns = target_rxns(~contains(target_rxns, '_REV'));

% find the genes in GPR rules of target rxns
enzyme_list = [];

for i = 1:numel(target_rxns)
    ind = find(ismember(auto_model.rxns, target_rxns{i}));
    t_enzyme = regexp(auto_model.rules{ind}, 'x\((\d+)\)', 'tokens');
    numbers = cellfun(@(x) str2double(x{1}), t_enzyme);
    tmp_cell = [];
    for j = 1:numel(numbers)
        tmp_cell = [tmp_cell, auto_model.genes(numbers(j))];
    end
    enzyme_list{i} = unique(tmp_cell);
end

enzyme_list = enzyme_list';

for i = 1:numel(enzyme_list)
    if iscell(enzyme_list{i})
        enzyme_list{i} = strjoin(enzyme_list{i}, ';');
    end
end

% find the rxns index
rxnIndex = rxn_cell(ismember(rxn_list, target_rxns));

% Create final table and write to file
final_table = cell2table([rxnIndex, target_rxns, enzyme_list], ...
    'VariableNames', {'RxnIndex', 'RxnID', 'Enzymes'});
writetable(final_table, 'Data/Reactions/list_of_rxns_with_at_least_one_enzyme_in_both.csv');

%% find the transport rxns and their associated enzyme
trans_index = find(contains(auto_model.rxnNames, 'transport'))

trans_index_with_GPR = [];

for i = 1:numel(trans_index)
    if ~isempty(auto_model.rules{trans_index(i)})
        if ~contains(auto_model.rxns(trans_index(i)), '_REV')
            trans_index_with_GPR = [trans_index_with_GPR; trans_index(i)];
        end
    end
end

transport_rxn = [];
associated_enzyme = [];

for i = 1:numel(trans_index_with_GPR)
    ind = trans_index_with_GPR(i);
    t_enzyme = regexp(auto_model.rules{ind}, 'x\((\d+)\)', 'tokens');
    numbers = cellfun(@(x) str2double(x{1}), t_enzyme);
    for j = 1:numel(numbers)
        transport_rxn = [transport_rxn; auto_model.rxns(ind)];
        associated_enzyme = [associated_enzyme; auto_model.genes(numbers(j))];
    end
end

final_table = cell2table([transport_rxn, associated_enzyme]);
final_table.Properties.VariableNames = {'RxnID', 'Enzymes'};
writetable(T, 'Data/Reactions/list_of_transport_rxn_with_enzyme.csv');