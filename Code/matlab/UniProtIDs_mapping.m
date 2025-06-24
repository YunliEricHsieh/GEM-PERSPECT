% translate the JGIv5.6 IDs to Uniport IDs
mutant_table = readtable('Data/Mutant_phenotypes_table_filtered_by_num_and_cor.csv');
table_for_map = readtable('Data/mart_Cre_Uniprot.txt');
 
% map IDs
for i = 1:numel(mutant_table.Gene)
    id = find(contains(table_for_map.GeneName, mutant_table.Gene{i}));

    if ~isempty(id)
        combinedString = [];

        for j = 1:numel(id)
            if j == 1
                combinedString = table_for_map.UniProtID{id(j)};
                mutant_table.UniProtID(i) = {combinedString};
            else 
                combinedString = [combinedString ' or ' table_for_map.UniProtID{id(j)}];
                mutant_table.UniProtID(i) = {combinedString};
            end
        end
    end
end

% remove the duplicated IDs
removeDuplicates = @(x) strjoin(unique(strsplit(x, ' or ')), ' or ');

for i = 1:numel(mutant_table.UniProtID)
    if ~isempty(mutant_table.UniProtID{i})
        mutant_table.UniProtID(i) = cellfun(removeDuplicates, mutant_table.UniProtID(i), 'UniformOutput', false);
    end
end

writetable(mutant_table, 'Data/Mutant_phenotypes_table_filtered_with_uniportIDs.csv')

