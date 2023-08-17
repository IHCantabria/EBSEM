function [S]= yates09(settings, params)

a = -params(1);
b = params(2);
cacr = -params(3);
cero = -params(4);


Seq = (settings.Y09.E - b) ./ a;

S = nan(size(settings.Y09.E));

S(1) = settings.Y09.S0;

for i = 1:numel(S)-1
    if S(i) < Seq(i+1)
        S(i+1) = ((S(i)-Seq(i+1)).*...
            exp(-1. * a *cacr *(settings.Y09.E(i+1) ^ 0.5)*settings.Y09.dt))+Seq(i+1);
    else
        S(i+1) = ((S(i)-Seq(i+1)).*...
            exp(-1. * a *cero *(settings.Y09.E(i+1) ^ 0.5)*settings.Y09.dt))+Seq(i+1);
    end
end
