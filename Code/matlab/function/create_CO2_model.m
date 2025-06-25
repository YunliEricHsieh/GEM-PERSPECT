function [CO2model, original_model] = create_CO2_model(model, alpha)

CO2model = model;
original_model = model;

% find the index for growth rate
bio_index = find(model.c == 1);

% find the optimal biomass value 
opt = optimizeCbModel(model, 'max');
opt_bio = opt.f;

% find index of uptake CO2 rxn
co2_index = find(ismember(model.rxns,'EX_co2_e_REV'));

% calculate the uptake flux value of optimal growth rate
% max flux 
model.lb(bio_index) = opt_bio;
model.ub(bio_index) = opt_bio;
opt1 = optimizeCbModel(model);
max_co2 = opt1.v(co2_index);

% min flux
model.lb(bio_index) = opt_bio*alpha;
model.ub(bio_index) = opt_bio*alpha;
opt2 = optimizeCbModel(model);
min_co2 = opt2.v(co2_index);

% fix the max CO2 for original model
original_model.ub(co2_index) = max_co2;

% increase 3% of CO2 for CO2 model
CO2model.ub(co2_index) = max_co2 + (max_co2-min_co2)*0.03;
