function [best_solution, min_value] = sce_ua2(f, x0, ngen, npop, npar, mag, lb, ub)
    % Initialize the population
    pop = zeros(npop, npar);
    pop(1, :) = x0;
    for i = 2:npop
        pop(i, :) = lb + (ub - lb) .* rand(1, npar);
    end

    fvals = zeros(1, npop);

    % Main loop
    for gen = 1:ngen

        % Evaluate the function for all individuals in the population
        parfor i = 1:npop
            fvals(i) = f(pop(i, :));
        end

        % Shuffle the population
        indices = randperm(npop);
        pop = pop(indices, :);
        fvals = fvals(indices);

        % Complex evolution
        for i = 1:2:npop
            y1 = fvals(i);
            y2 = fvals(i + 1);

            if y1 < y2
                x1 = pop(i, :);
                x2 = pop(i + 1, :);
                z = x1 + 0.5 .* (x2 - x1) + 0.5 .* mag .* (2. * rand(1, npar) - 1);
                z = max(z, lb); % clip to lower bounds
                z = min(z, ub); % clip to upper bounds
                pop(i + 1, :) = z;
            else
                x2 = pop(i, :);
                x1 = pop(i + 1, :);
                z = x2 + 0.5 .* (x1 - x2) + 0.5 .* mag .* (2. * rand(1, npar) - 1);
                z = max(z, lb); % clip to lower bounds
                z = min(z, ub); % clip to upper bounds
                pop(i, :) = z;
            end
        end

        if mod(gen, floor(ngen / 20)) == 0
            fprintf('Progress = %.2f %%\n', 100 * gen / ngen);
            fprintf('Generation: %d / %d\n', gen, ngen);
        end

        % Check convergence
        if gen > 5
            if all(abs(fvals(1) - fvals) < 1e-4)
                fprintf('Converged at generation %d.\n', gen);
                break;
            end
        end

    end

    % Return the best solution
    [~, best_index] = min(fvals);
    best_solution = pop(best_index, :);
    min_value = fvals(best_index);
end
