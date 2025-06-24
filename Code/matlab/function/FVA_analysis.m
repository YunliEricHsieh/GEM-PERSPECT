function [minFlux, maxFlux] = FVA_analysis (model, percentage, rxn_list)

minFlux = cell(numel(rxn_list), 1);
maxFlux = cell(numel(rxn_list), 1);

bio_index = find(model.c == 1);

% find the optimal biomass value 
opt = optimizeCbModel(model, 'max');
opt_bio = opt.f;

model.lb(bio_index) = opt_bio*(percentage/100);
model.ub(bio_index) = opt_bio*(percentage/100);

problem.A = sparse(model.S); 
problem.rhs = model.b;
problem.lb = model.lb;
problem.ub = model.ub;
problem.sense = repelem('=',size(model.b,1),1);
problem.vtype = repelem('C',size(model.S,2),1);
problem.obj = zeros(size(model.S,2), 1);
params.Threads = 1;

parfor i = 1:numel(rxn_list)

    model_local = model;
    problem_local = problem;

    rxn_index = find(ismember(model_local.rxns, rxn_list{i}));
    problem_local.obj(rxn_index) = 1;

    % Solve using Gurobi - 'max'
    problem_local.modelsense = 'max';
    flux_solution = gurobi(problem_local, params);

    if isfield(flux_solution, 'objval')
        maxFlux{i} = flux_solution.objval;
    else
        maxFlux{i} = 'Inf';
    end

    % Solve using Gurobi - 'min'
    problem_local.modelsense = 'min';
    flux_solution = gurobi(problem_local, params);

    if isfield(flux_solution, 'objval')
        minFlux{i} = flux_solution.objval;
    else
        minFlux{i} = 'Inf';
    end     

end

