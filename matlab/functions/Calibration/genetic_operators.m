function offspring = genetic_operators(parents)
[num_offspring, npar] = size(parents);
offspring = zeros(num_offspring, npar);

for i = 1:2:num_offspring
    % Select parents to breed
    p1 = parents(i, :);
    p2 = parents(i + 1, :);
    
    % Apply crossover operator
    [o1, o2] = crossover(p1, p2);
    
    % Apply mutation operator
    o1 = mutation(o1);
    o2 = mutation(o2);
    
    offspring(i, :) = o1;
    offspring(i + 1, :) = o2;
end

