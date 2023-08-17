function y = JV_APPROX(v, F, x)

% calculate perturbation
dim = size(x, 2);


sum=cumsum(sqrt(eps) .* (1 + x));

if norm(v, 2) > eps
    
    per = (1 ./ (dim .* norm(v, 2))) .* sum;
    
else

    per = sum ./ dim;
    
end

R = F(x); % unperturbed residual

xper = x + per .* v; % perturbed vector

Rper = F(xper); % perturbed residual

y = (Rper - R) ./ per; % approximation of jacobian action on krylov vector

if ~iscolumn(y)
    y=y';
end

end