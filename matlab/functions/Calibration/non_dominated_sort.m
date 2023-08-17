function [NEW_pop, NEW_obj] = non_dominated_sort(population, objectives)
% Sort population into non-dominated fronts
fronts = sort_fronts(objectives);

% Select individuals from each front
selected = select_individuals(fronts);

% Update population and objectives arrays
NEW_pop = population(selected, :);
NEW_obj = objectives(selected, :);
