function selected = select_individuals(fronts)
% Select individuals from each front
selected = [];
for i = 1:length(fronts)
    front = fronts{i};
    n = length(front);
    if n > 0
        if n <= length(selected)
            selected(1:n) = front;
        else
            selected = [selected; front'];
        end
    end
end
   

