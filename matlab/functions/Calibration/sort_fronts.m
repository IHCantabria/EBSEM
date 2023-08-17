function fronts = sort_fronts(objectives)
% Sort population into non-dominated fronts
npop = size(objectives, 1);

dominated_by = zeros(npop, 1);
n_dominated = zeros(npop, 1);

fronts = cell(1, 1);
fronts{1} = [];

for p = 1:npop
    for q = 1:npop
        if p == q
            continue;
        end

        if dominates(objectives(p, :), objectives(q, :))
            dominated_by(q) = dominated_by(q) + 1;
        elseif dominates(objectives(q, :), objectives(p, :))
            n_dominated(p) = n_dominated(p) + 1;
        end
    end

    if n_dominated(p) == 0
        fronts{1} = [fronts{1}, p];
    end
end

i_front = 1;
while ~isempty(fronts{i_front})
    next_front = [];
    for p = fronts{i_front}
        for q = find(dominated_by == p)'
            n_dominated(q) = n_dominated(q) - 1;
            if n_dominated(q) == 0
                next_front = [next_front, q];
            end
        end
    end
    fronts{i_front + 1} = next_front;
    i_front = i_front + 1;
end

% Remove the last empty front
fronts = fronts(1:end-1);


