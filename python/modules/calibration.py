import numpy as np

def sce_ua2(f, x0, ngen, npop, npar, mag, lb, ub):
    # Initialize the population
    pop = np.zeros((npop, npar))
    pop[0, :] = x0
    for i in range(1, npop):
        pop[i, :] = lb + (ub - lb) * np.random.rand(1, npar)

    fvals = np.zeros(npop)
    
    # Main loop
    for gen in range(1, ngen+1):

        # Evaluate the function for all individuals in the population
        func = np.vectorize(lambda i: f(pop[i,:]))
        fvals = func(np.arange(0,npop),)
        # for i in range(npop):
        #     fvals[i] = f(pop[i, :])
        
        # Shuffle the population
        indices = np.random.permutation(npop)
        pop = pop[indices, :]
        fvals = fvals[indices]

        # Complex evolution
        for i in range(0, npop, 2):
            y1 = fvals[i]
            y2 = fvals[i + 1]

            if y1 < y2:
                x1 = pop[i, :]
                x2 = pop[i + 1, :]
                z = x1 + 0.5 * (x2 - x1) + 0.5 * mag * (2. * np.random.rand(1, npar) - 1)
                z = np.maximum(z, lb)  # clip to lower bounds
                z = np.minimum(z, ub)  # clip to upper bounds
                pop[i + 1, :] = z
            else:
                x2 = pop[i, :]
                x1 = pop[i + 1, :]
                z = x2 + 0.5 * (x1 - x2) + 0.5 * mag * (2. * np.random.rand(1, npar) - 1)
                z = np.maximum(z, lb)  # clip to lower bounds
                z = np.minimum(z, ub)  # clip to upper bounds
                pop[i, :] = z

        if gen % int(ngen / 20) == 0:
            print(f'Progress = {100 * gen / ngen:.2f} %')
            print(f'Generation: {gen} / {ngen}')
            print(f'Metric Value: {np.min(fvals):.2f}')
        # Check convergence
        if gen > 5:
            if np.all(np.abs(fvals[0] - fvals) < 1e-4):
                print(f'Converged at generation {gen}.')
                break

    # Return the best solution
    best_index = np.argmin(fvals)
    best_solution = pop[best_index, :]
    min_value = fvals[best_index]
    
    return best_solution, min_value

