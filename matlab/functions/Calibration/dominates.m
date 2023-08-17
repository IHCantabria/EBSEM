function dom = dominates(obj1, obj2)
% Determine if obj1 dominates obj2
% A solution dominates another if it is better in all objectives
% and not worse in at least one objective.

dominates = true;
worse = false;

for i = 1:length(obj1)
    if obj1(i) > obj2(i)
        dominates = false;
        break;
    elseif obj1(i) < obj2(i)
        worse = true;
    end
end

dom = dominates && worse;