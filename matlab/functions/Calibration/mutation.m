function mutated = mutation(o)
% Apply mutation operator to offspring
% Add Gaussian noise with mean 0 and standard deviation 1 to each gene
mutated = o + randn(size(o)) .* 0.2;
