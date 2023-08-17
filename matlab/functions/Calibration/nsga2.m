function [population, objectives] = nsga2(obj_fun, npar, ngen, npop, ub, lb, nobj)
% Initialize population of candidate solutions
population = randn(npop, npar);

objectives = zeros(npop, nobj);
% Evaluate initial population
parfor i = 1:npop
    scaled_individual = ((population(i, :) - (-1)) ./ (1 - (-1))) .* (ub - lb) + lb;
    [a, b, c] = obj_fun(scaled_individual);
    objectives(i, :) = [a, b, c]; 
end

% Initialize iteration counter
t = 1;

while t <= ngen
    % Select parents for breeding
    parents = tournament_selection(population, objectives);
    
    % Create offspring using genetic operators
    offspring = genetic_operators(parents);
    noff = size(offspring, 1);
    
    offspring_obj = zeros(noff, nobj);
    % Evaluate offspring
    parfor i = 1:noff
        scaled_offspring = ((offspring(i, :) - (-1)) ./ (1 - (-1))) .* (ub - lb) + lb;
        [a, b, c] = obj_fun(scaled_offspring);
        offspring_obj(i, :) = [a, b, c];
    end
    
    % Update population with offspring and original solutions
    [population, objectives] = non_dominated_sort(...
        [population; offspring], ...
        [objectives; offspring_obj] ...
    );
    % Check length of sorted population and add random individuals if necessary
    pop_len = size(population, 1);
    if pop_len < npop && t ~= ngen
        num_missing = npop - pop_len;
        new_pop = randn(num_missing, npar);
        new_obj = zeros(num_missing, nobj);
        parfor i = 1:num_missing
            scaled_new_individual = ((new_pop(i, :) - (-1)) ./ (1 - (-1))) .* (ub - lb) + lb;
            [a, b, c] = obj_fun(scaled_new_individual);
            new_obj(i, :) = [a, b, c];
        end
        population = [population; new_pop];
        objectives = [objectives; new_obj];
    elseif pop_len > npop
        population = population(1:npop, :);
        objectives = objectives(1:npop, :);
    end
    
    % Update iteration counter
    if mod(t, floor(ngen / 20)) == 0
        fprintf('Progress = %.2f %%\n', 100 * t / ngen);
        fprintf('Generation: %d / %d\n', t, ngen);
    end
    t = t + 1;
end

for i = 1:size(population, 1)
    scaled_individual = ((population(i, :) - (-1)) ./ (1 - (-1))) .* (ub - lb) + lb;
    population(i, :) = scaled_individual;
end
