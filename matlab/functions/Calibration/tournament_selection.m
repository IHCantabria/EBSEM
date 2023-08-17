function selected = tournament_selection(population, objectives)
[npop, npar] = size(population);
num_selected = floor(npop / 2);
selected = zeros(num_selected, npar);

for i = 1:num_selected
    competitors = randperm(npop, 2);
    if dominates(objectives(competitors(1), :), objectives(competitors(2), :))
        selected(i, :) = population(competitors(1), :);
    else
        selected(i, :) = population(competitors(2), :);
    end
end