function [o1, o2] = crossover(p1, p2)
alpha = rand();
o1 = alpha * p1 + (1 - alpha) * p2;
o2 = (1 - alpha) * p1 + alpha * p2;

