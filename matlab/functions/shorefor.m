function [S]=shorefor(settings, params)

phi = params(1);
c = params(2);

D = 2*phi;

ii = (0:settings.SF.dt:(D-1)*24)';

phivecP = 10 .^ (-abs(ii) ./ (phi * 24));

IDX = length(phivecP);

phivecP = [zeros(IDX,1); phivecP];

vent = phivecP ./ sum(phivecP);

OmegaEQ = imfilter(settings.Omega - mean(settings.Omega), vent, 'same')...
    + mean(settings.Omega);

F = (settings.SF.P.^0.5) .* (OmegaEQ - settings.Omega) ./ std(OmegaEQ);


F(1:IDX-1) = 0;

S = nan(length(settings.Omega),1);
rero = F < 0;
racr = F >= 0;
S(1) = params(3);

r = abs(sum(F(racr)) ./ sum(F(rero)));

r_rero_F = r .* rero(2:end) .* F(2:end);
racr_F = racr(2:end) .* F(2:end);
r_rero_F_prev = r .* rero(1:end-1) .* F(1:end-1);
racr_F_prev = racr(1:end-1) .* F(1:end-1);
S(2:end) = 0.5 .* settings.SF.dt .* c .* cumsum(r_rero_F + racr_F + r_rero_F_prev...
    + racr_F_prev) + S(1);
