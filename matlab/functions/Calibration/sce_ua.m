function [best_solution, min_value] = sce_ua(f, x0, ngen, npop, npar, mag)
    % Initialize the population
    pop = zeros(npop, npar);
    for i = 1:npop
        pop(i, :) = x0 + mag .* (2. * rand(1, npar) - 1);
        % for j = 1:npar
        %     pop(i, j) = min(upper_bounds(j), max(lower_bounds(j), pop(i, j)));
        % end
    end

    % Main loop
    fvals = zeros(1, npop);
    for gen = 1:ngen
        % Shuffle the population
        for i = 1:npar
            pop(:, i) = pop(randperm(npop), i);
        end

        % Evaluate the function for all individuals in the population
        parfor i = 1:npop
            fvals(i) = f(pop(i, :));
        end

        % Complex evolution
        for i = 1:2:npop
            x1 = pop(i, :);
            x2 = pop(i + 1, :);
            y1 = fvals(i);
            y2 = fvals(i + 1);

            if y1 < y2
                % Use x1 as the better individual
                pop(i + 1, :) = x1 + 0.5 .* (x2 - x1) + 0.5 .* mag .* (2. * rand(1, npar) - 1);
            else
                % Use x2 as the better individual
                pop(i, :) = x2 + 0.5 .* (x1 - x2) + 0.5 .* mag .* (2. * rand(1, npar) - 1);
            end
        end
        if ceil(gen/100)-gen/100 == 0
            fprintf('Progress = %.2f %%\n', gen / ngen * 100);
            fprintf('Generation: %d / %d\n\n', gen, ngen);
        end
    end

    % Return the best solution
    [~, best_index] = min(fvals);
    best_solution = pop(best_index, :);
    min_value = fvals(best_index);
end