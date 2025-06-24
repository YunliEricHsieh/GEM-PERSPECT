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

auto_opt = optimizeCbModel(auto_model);
hetero_opt = optimizeCbModel(hetero_model);
mixo_opt = optimizeCbModel(mixo_model);

auto_index = find(auto_model.c == 1);
hetero_index = find(hetero_model.c == 1);
mixo_index = find(mixo_model.c == 1);

auto_model.lb(auto_index) = 0.9*auto_opt.f;
hetero_model.lb(hetero_index) = 0.9*hetero_opt.f;
mixo_model.lb(mixo_index) = 0.9*mixo_opt.f;

auto_model.ub(auto_index) = 0.9*auto_opt.f;
hetero_model.ub(hetero_index) = 0.9*hetero_opt.f;
mixo_model.ub(mixo_index) = 0.9*mixo_opt.f;

[a_sampleStruct, a_mixedFraction] = gpSampler(auto_model, 5000);
[h_sampleStruct, h_mixedFraction] = gpSampler(hetero_model, 5000);
[m_sampleStruct, m_mixedFraction] = gpSampler(mixo_model, 5000);

%% create the rxn index
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
re_rxn_cells = [forward_cells; rerxn_cell];
re_rxn_list = [forward_rxns; rerxn];

% filter out the reversible reactions
rxn_list = rxn_list(~ismember(rxn_list, re_rxn_list));
rxn_cell = rxn_cell(~ismember(rxn_cell, re_rxn_cells));

%% find the flux distribution
% auto_model
meanFlux = cell(height(rxn_cell),1);
stdFlux = cell(height(rxn_cell),1);

a_fluxSamples = a_sampleStruct.points;
for i = 1:height(rxn_cell)
    rxnIndex = find(ismember(auto_model.rxns, rxn_list{i}));

    targetFluxDistribution = a_fluxSamples(rxnIndex, :);
    meanFlux{i} = mean(targetFluxDistribution);
    stdFlux{i} = std(targetFluxDistribution);

    clear targetFluxDistribution
end

tmp_table = cell2table([rxn_cell, rxn_list, meanFlux, stdFlux]);
tmp_table.Properties.VariableNames = {'RxnIndex', 'RxnID', 'meanFlux', 'stdFlux'};

writetable(tmp_table, 'Results/flux_sampling/auto_sampling.csv');

% hetero_model
meanFlux = cell(height(rxn_cell),1);
stdFlux = cell(height(rxn_cell),1);

h_fluxSamples = h_sampleStruct.points;
for i = 1:height(rxn_cell)

    rxnIndex = find(ismember(hetero_model.rxns, rxn_list{i}));

    targetFluxDistribution = h_fluxSamples(rxnIndex, :);
    meanFlux{i} = mean(targetFluxDistribution);
    stdFlux{i} = std(targetFluxDistribution);

    clear targetFluxDistribution
end

tmp_table = cell2table([rxn_cell, rxn_list, meanFlux, stdFlux]);
tmp_table.Properties.VariableNames = {'RxnIndex', 'RxnID', 'meanFlux', 'stdFlux'};

writetable(tmp_table, 'Results/flux_sampling/hetero_sampling.csv');

% mixo_model
meanFlux = cell(height(rxn_cell),1);
stdFlux = cell(height(rxn_cell),1);

m_fluxSamples = m_sampleStruct.points;
for i = 1:height(rxn_cell)

    rxnIndex = find(ismember(mixo_model.rxns, rxn_list{i}));

    targetFluxDistribution = m_fluxSamples(rxnIndex, :);
    meanFlux{i} = mean(targetFluxDistribution);
    stdFlux{i} = std(targetFluxDistribution);

    clear targetFluxDistribution
end

tmp_table = cell2table([rxn_cell, rxn_list, meanFlux, stdFlux]);
tmp_table.Properties.VariableNames = {'RxnIndex', 'RxnID', 'meanFlux', 'stdFlux'};

writetable(tmp_table, 'Results/flux_sampling/mixo_sampling.csv');

%% find the flux distribution for reversible reactions
% auto_model
meanFlux = cell(height(rerxn),1);
stdFlux = cell(height(rerxn),1);

a_fluxSamples = a_sampleStruct.points;
for i = 1:height(rerxn)
    rxnIndex_b = find(ismember(auto_model.rxns, rerxn{i}));
    rxnIndex_f = find(ismember(auto_model.rxns, strrep(rerxn{i}, '_REV', '')));

    targetFluxDistribution_f = a_fluxSamples(rxnIndex_f, :);
    targetFluxDistribution_b = a_fluxSamples(rxnIndex_b, :);

    meanFlux{i} = mean(abs(targetFluxDistribution_f - targetFluxDistribution_b));
    stdFlux{i} = std(abs(targetFluxDistribution_f - targetFluxDistribution_b));

    clear targetFluxDistribution_f targetFluxDistribution_b
end

tmp_table = cell2table([forward_cells, forward_rxns, meanFlux, stdFlux]);
tmp_table.Properties.VariableNames = {'RxnIndex', 'RxnID', 'meanFlux', 'stdFlux'};

writetable(tmp_table, 'Results/flux_sampling/auto_sampling_re.csv');

% hetero_model
meanFlux = cell(height(rerxn),1);
stdFlux = cell(height(rerxn),1);

h_fluxSamples = h_sampleStruct.points;
for i = 1:height(rerxn)
    rxnIndex_b = find(ismember(hetero_model.rxns, rerxn{i}));
    rxnIndex_f = find(ismember(hetero_model.rxns, strrep(rerxn{i}, '_REV', '')));

    targetFluxDistribution_f = h_fluxSamples(rxnIndex_f, :);
    targetFluxDistribution_b = h_fluxSamples(rxnIndex_b, :);

    meanFlux{i} = mean(abs(targetFluxDistribution_f - targetFluxDistribution_b));
    stdFlux{i} = std(abs(targetFluxDistribution_f - targetFluxDistribution_b));

    clear targetFluxDistribution_f targetFluxDistribution_b
end

tmp_table = cell2table([forward_cells, forward_rxns, meanFlux, stdFlux]);
tmp_table.Properties.VariableNames = {'RxnIndex', 'RxnID', 'meanFlux', 'stdFlux'};

writetable(tmp_table, 'Results/flux_sampling/hetero_sampling_re.csv');

% mixo_model
meanFlux = cell(height(rerxn),1);
stdFlux = cell(height(rerxn),1);

m_fluxSamples = m_sampleStruct.points;
for i = 1:height(rerxn)
    rxnIndex_b = find(ismember(mixo_model.rxns, rerxn{i}));
    rxnIndex_f = find(ismember(mixo_model.rxns, strrep(rerxn{i}, '_REV', '')));

    targetFluxDistribution_f = m_fluxSamples(rxnIndex_f, :);
    targetFluxDistribution_b = m_fluxSamples(rxnIndex_b, :);

    meanFlux{i} = mean(abs(targetFluxDistribution_f - targetFluxDistribution_b));
    stdFlux{i} = std(abs(targetFluxDistribution_f - targetFluxDistribution_b));

    clear targetFluxDistribution_f targetFluxDistribution_b
end

tmp_table = cell2table([forward_cells, forward_rxns, meanFlux, stdFlux]);
tmp_table.Properties.VariableNames = {'RxnIndex', 'RxnID', 'meanFlux', 'stdFlux'};

writetable(tmp_table, 'Results/flux_sampling/mixo_sampling_re.csv');