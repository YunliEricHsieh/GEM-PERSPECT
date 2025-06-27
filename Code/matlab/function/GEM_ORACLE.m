function flux_values = GEM_ORACLE(models, ratios, rxn_list, alpha, type)

indices = ~isnan(ratios);
model_index = zeros(8, 1);

if indices(1) == 1
    model_index(1) = 1;
    model_index(2) = 1;
end

if indices(2) == 1
    model_index(1) = 1;
    model_index(3) = 1;
end

if indices(3) == 1
    model_index(2) = 1;
    model_index(3) = 1;
end

if indices(4) == 1
    model_index(1) = 1;
    model_index(4) = 1;
end

if indices(5) == 1
    model_index(2) = 1;
    model_index(5) = 1;
end

if indices(6) == 1
    model_index(2) = 1;
    model_index(6) = 1;
end

if indices(7) == 1
    model_index(2) = 1;
    model_index(7) = 1;
end

if indices(8) == 1
    model_index(2) = 1;
    model_index(8) = 1;
end

% find the testing models
test_models = models(logical(model_index));

%% prepare the matrix
numModels = numel(test_models);
bio_index = NaN(8,1);
opt_bio = NaN(8,1);
ub = [];
lb = [];
nrxn = zeros(8,1);
nmet = zeros(8,1);

ind = find(model_index);
for i = 1:numModels
    % find the index of biomass rxn
    bio_index(ind(i)) = find(test_models{i}.c == 1);

    % find the optimal biomass value 
    opt = optimizeCbModel(test_models{i}, 'max');
    opt_bio(ind(i)) = opt.f;

    % growth rate constraint
    test_models{i}.lb(bio_index(ind(i))) = alpha*opt_bio(ind(i));

    % create ub and lb 
    ub = [ub; test_models{i}.ub];
    lb = [lb; test_models{i}.lb];

    % calculate the total rxn and met number
    nrxn(ind(i)) = numel(test_models{i}.rxns);
    nmet(ind(i)) = numel(test_models{i}.mets);
end

% build the matrix
% Preallocate Aeq matrix
Aeq = zeros(sum(nmet), sum(nrxn)); 
Aeq1 = zeros(numModels-1, size(Aeq,2));  

row = 0;
col = 0;

for i = 1:numModels
    Aeq(row + (1:nmet(ind(i))), col + (1:nrxn(ind(i)))) = test_models{i}.S;
    row = row + nmet(ind(i));
    col = col + nrxn(ind(i));
end

beq = zeros(size(Aeq,1), 1);

% v_a - v_m = 0
beq1 = zeros(numModels-1,1);

% v_bio - v_bio*r < w and -v_bio + v_bio*r < w (w = 0.05*r)
w = zeros(8, 1);
Aineq = zeros(16, size(Aeq,2));

if indices(1)
    Aineq(1, [bio_index(1), nrxn(1)+bio_index(2)]) = [1, -ratios(1)];
    Aineq(9, [bio_index(1), nrxn(1)+bio_index(2)]) = [-1, ratios(1)];
    w(1) = ratios(1)*0.05;
end

if indices(2)
    Aineq(2, [bio_index(1), nrxn(1)+nrxn(2)+bio_index(3)]) = [1, -ratios(2)];
    Aineq(10, [bio_index(1), nrxn(1)+nrxn(2)+bio_index(3)]) = [-1, ratios(2)];
    w(2) = ratios(2)*0.05;
end

if indices(3)
    Aineq(3, [nrxn(1)+bio_index(2), nrxn(1)+nrxn(2)+bio_index(3)]) = [1, -ratios(3)];
    Aineq(11, [nrxn(1)+bio_index(2), nrxn(1)+nrxn(2)+bio_index(3)]) = [-1, ratios(3)];
    w(3) = ratios(3)*0.05;
end

if indices(4)
    Aineq(4, [bio_index(1), nrxn(1)+nrxn(2)+nrxn(3)+bio_index(4)]) = [1, -ratios(4)];
    Aineq(12, [bio_index(1), nrxn(1)+nrxn(2)+nrxn(3)+bio_index(4)]) = [-1, ratios(4)];
    w(4) = ratios(4)*0.05;
end

if indices(5)
    Aineq(5, [nrxn(1)+bio_index(2), nrxn(1)+nrxn(2)+nrxn(3)+nrxn(4)+bio_index(5)]) = [1, -ratios(5)];
    Aineq(13, [nrxn(1)+bio_index(2), nrxn(1)+nrxn(2)+nrxn(3)+nrxn(4)+bio_index(5)]) = [-1, ratios(5)];
    w(5) = ratios(5)*0.05;
end

if indices(6)
    Aineq(6, [nrxn(1)+bio_index(2), nrxn(1)+nrxn(2)+nrxn(3)+nrxn(4)+nrxn(5)+bio_index(6)]) = [1, -ratios(6)];
    Aineq(14, [nrxn(1)+bio_index(2), nrxn(1)+nrxn(2)+nrxn(3)+nrxn(4)+nrxn(5)+bio_index(6)]) = [-1, ratios(6)];
    w(6) = ratios(6)*0.05;
end

if indices(7)
    Aineq(7, [nrxn(1)+bio_index(2), nrxn(1)+nrxn(2)+nrxn(3)+nrxn(4)+nrxn(5)+nrxn(6)+bio_index(7)]) = [1, -ratios(7)];
    Aineq(15, [nrxn(1)+bio_index(2), nrxn(1)+nrxn(2)+nrxn(3)+nrxn(4)+nrxn(5)+nrxn(6)+bio_index(7)]) = [-1, ratios(7)];
    w(7) = ratios(7)*0.05;
end

if indices(8)
    Aineq(8, [nrxn(1)+bio_index(2), nrxn(1)+nrxn(2)+nrxn(3)+nrxn(4)+nrxn(5)+nrxn(6)+nrxn(7)+bio_index(8)]) = [1, -ratios(8)];
    Aineq(16, [nrxn(1)+bio_index(2), nrxn(1)+nrxn(2)+nrxn(3)+nrxn(4)+nrxn(5)+nrxn(6)+nrxn(7)+bio_index(8)]) = [-1, ratios(8)];
    w(8) = ratios(8)*0.05;
end

% remove the empty rows
Aineq(~any(Aineq, 2), :) = [];

bineq = [w(indices); w(indices)];

% build the matrix for gurobi solver
problem.A = sparse([Aeq; Aineq; Aeq1]); 
problem.rhs = [beq; bineq; beq1];
problem.lb = lb;
problem.ub = ub;
problem.vtype = repelem('C',size(Aeq,2),1);
problem.sense = [repelem('=',size(beq,1),1); repelem('<',size(bineq,1),1); repelem('=',size(beq1,1),1)];

if ismember(type, {'max'})
    problem.modelsense = 'max';
elseif ismember(type, {'min'})
    problem.modelsense = 'min';
end        

%% loop all the reactions
params.Threads = 2;

flux_values = cell(1, numel(rxn_list));

parfor j = 1:numel(rxn_list) 

    % Local copies of variables inside parfor loop
    problem_local = problem;
      
    % find the target rxns
    target_rxns = zeros(numel(test_models), 1);

    for i = 1:numel(test_models)
        target_rxns(i) = find(ismember(test_models{i}.rxns, rxn_list{j}));
    end

    % assign the reaction constraint into matrix
    % v_a - v_m = 0
    col = 0;
    rowT = size(problem.A,1);

    for i = 1:numModels-1
        row = rowT-(numModels-1)+i;
        problem_local.A(row, [col+target_rxns(i), col+nrxn(ind(i))+target_rxns(i+1)]) = [1, -1];
        col = col + nrxn(ind(i));
    end

    % calculate the min flux sum value across all the reactions
    problem1 = problem_local;
    problem1.modelsense = 'min';
    problem1.obj = ones(size(Aeq, 2), 1);

    % Solve min flux sum 
    min_flux_sum = gurobi(problem1, params);

    if isfield(min_flux_sum, 'objval')
        % add the min flux sum as constraint
        problem_local.A = sparse([problem_local.A; ones(1, size(Aeq, 2))]);
        problem_local.rhs = [problem_local.rhs; min_flux_sum.objval];
        problem_local.sense = [problem_local.sense; '='];

        % determine the objective 
        obj = zeros(size(Aeq,2), 1);
        obj(target_rxns(1)) = 1;

        problem_local.obj = obj;

        % Solve using Gurobi
        flux_solution = gurobi(problem_local, params);

        if isfield(flux_solution, 'objval')
            flux_values{j} = flux_solution.objval;
        else
            flux_values{j} = 'Inf';
        end
    else
        flux_values{j} = 'Inf';
    end

end























