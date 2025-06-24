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

%% run FVA
ncpu = 20;
parpool(ncpu);
parfevalOnAll(@maxNumCompThreads, 0, 1); 

% FVA analysis
[auto_minFlux, auto_maxFlux] = FVA_analysis (auto_model, 10, rxn_list);
[mixo_minFlux, mixo_maxFlux]  = FVA_analysis (mixo_model, 10, rxn_list);
[hetero_minFlux, hetero_maxFlux] = FVA_analysis (hetero_model, 10, rxn_list);

% save tables
auto_FVA = table(rxn_list, auto_minFlux, auto_maxFlux, 'VariableNames', {'RxnID', 'minFlux', 'maxFlux'});
mixo_FVA = table(rxn_list, mixo_minFlux, mixo_maxFlux, 'VariableNames', {'RxnID', 'minFlux', 'maxFlux'});
hetero_FVA = table(rxn_list, hetero_minFlux, hetero_maxFlux, 'VariableNames', {'RxnID', 'minFlux', 'maxFlux'});

writetable(auto_FVA,'Results/FVA/auto_FVA_10p.csv');
writetable(mixo_FVA,'Results/FVA/mixo_FVA_10p.csv');
writetable(hetero_FVA,'Results/FVA/hetero_FVA_10p.csv');