%% load model
model = readCbModel('Data/pciCre1355/NDLadpraw_Autotrophic_Rep1.xml');

% find the associated substrates for the target reaction
rxn_list = {'ATAH','LLDHm','ATAH_co2','PLDAGAT1819Z18111Z18111Z1','URCB'};

rxnIndices    = findRxnIDs(model, rxn_list);
isReversible = ismember(strcat(rxn_list, '_REV'), model.rxns);

RxnID       = {};
BiGG_ID     = {};
MetFormula  = {};
Substrate   = {};
Pathway     = {};

pattern = '->'; 
coefRE  = '^\d+\s+';

for i = 1:numel(rxn_list)
    rxnID    = rxn_list{i};
    ridx     = rxnIndices(i);
    pathway  = model.subSystems{ridx};
    formula  = printRxnFormula(model, rxnID);
    sides    = strtrim(strsplit(formula{1}, pattern));
    LHSmets  = regexprep(strtrim(strsplit(sides{1}, '+')), coefRE, '');
    allSides = LHSmets(:);
    
    % include RHS if reversible
    if isReversible(i)
        RHSmets = regexprep(strtrim(strsplit(sides{2}, '+')), coefRE, '');
        allSides = [allSides; RHSmets(:)];
    end
    
    % lookup metabolite indices
    metIdx = cellfun(@(m)find(strcmp(model.mets, m),1), allSides);

    % append data
    n = numel(allSides);
    RxnID      = [RxnID;       repmat({rxnID}, n,1)];
    BiGG_ID    = [BiGG_ID;     allSides];
    MetFormula = [MetFormula;  model.metFormulas(metIdx)];
    Substrate  = [Substrate;   model.metNames(metIdx)];
    Pathway    = [Pathway;     repmat({pathway}, n,1)];
end

% Build and write final table
T = table(RxnID, BiGG_ID, Substrate, MetFormula, Pathway, ...
    'VariableNames',{'RxnID','BiGG_ID','Substrate','MetFormula','Pathway'});
writetable(T,'Results/CataPro/enzymes_and_substrates.csv');

%% find the reaction in the target pathway
targetsubSystem = {'Purine metabolism', 'Pyruvate metabolism', 'Glycerolipid metabolism', 'Urea degradation'};

RxnID       = {};
BiGG_ID     = {};
Enzyme   = {};
MetFormula  = {};
Substrate   = {};
Pathway     = {};

pattern1 = 'x\((\d+)\)';

for i = 1:numel(targetsubSystem)
    % find the rxn indexes for the target subsystem
    index = [];
    for j = 1:numel(model.subSystems)
        if ~isempty(find(contains(model.subSystems{j}, targetsubSystem{i})))
            index = [index; j];
        end
    end
    
    rxn_list = model.rxns(index);
    % remove sudo reactions
    rxn_list   = rxn_list(~contains(rxn_list, 'No'));
    rxnIndices = findRxnIDs(model, rxn_list);

    for j = 1:numel(rxn_list)
        ridx   = rxnIndices(j);
        % find the rxn wiht GPR rules
        if ~isempty(model.rules{rxnIndices(j)})
            rxnID      = rxn_list{j};
            formula    = printRxnFormula(model, rxnID);
            sides      = strtrim(strsplit(formula{1}, pattern));
            LHSmets    = regexprep(strtrim(strsplit(sides{1}, '+')), coefRE, '');
            allSides   = LHSmets(:);
            % lookup metabolite indices
            metIdx     = cellfun(@(m)find(strcmp(model.mets, m),1), allSides);
            % find associated enzyme
            eidx       = regexp(model.rules{rxnIndices(j)}, pattern1, 'tokens');
            numbers    = cellfun(@(x) str2double(x{1}), eidx)
            enzyme = unique(model.genes(numbers));
            % append data
            n          = numel(allSides);
            m          = numel(enzyme);
            RxnID      = [RxnID;       repmat({rxnID}, n*m,1)];
            BiGG_ID    = [BiGG_ID;     repmat(allSides, m, 1)];
            MetFormula = [MetFormula;  repmat(model.metFormulas(metIdx), m, 1)];
            Substrate  = [Substrate;   repmat(model.metNames(metIdx), m, 1)];
            Pathway    = [Pathway;     repmat({targetsubSystem{i}}, n*m,1)];
            Enzyme     = [Enzyme;      repelem(enzyme, n,1)];
        end
    end
end

% Build and write final table
T = table(RxnID, BiGG_ID, Substrate, MetFormula, Pathway, Enzyme, ...
    'VariableNames',{'RxnID','BiGG_ID','Substrate','MetFormula','Pathway', 'Enzyme'});
writetable(T,'Results/CataPro/enzymes_and_substrates_for_reference.csv');