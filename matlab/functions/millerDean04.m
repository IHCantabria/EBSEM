function [Y]= millerDean04(settings, params)

DY0 = params(1);
cacr = params(2);
cero = params(3);

kero = zeros(length(settings.Hb),1);
kacr = zeros(length(settings.Hb),1);

if settings.MD.flagP == 1
    kero(:,:) = cero;
    kacr(:,:) = cacr;
elseif settings.MD.flagP == 2
    kero(:,:) = cero.* settings.Hb .^ 2;
    kacr(:,:) = cacr.* settings.Hb .^ 2;
elseif settings.MD.flagP == 3
    kero(:,:) = cero.* settings.Hb .^ 3;
    kacr(:,:) = cacr.* settings.Hb .^ 3;
elseif settings.MD.flagP == 4
    kero(:,:) = cero.* settings.Omega;
    kacr(:,:) = cacr.* settings.Omega;
end

Y = nan(length(settings.Hb),1);
wl = 0.106 .* settings.Hb + settings.MD.sl;
Wast = wast(settings.depthb, settings.d50);
yeq = DY0 - Wast .* wl ./ (settings.MD.Hberm + settings.depthb);
Y(1)=settings.MD.Yi;

for i = 2:length(settings.MD.sl)
    r = yeq(i) - Y(i-1) > 0;
    k = kacr(i) * r + kero(i) * ~r;
    A = k .* settings.MD.dt .* 0.5;
    Y(i) = (Y(i-1) + A .* (yeq(i) + yeq(i-1) - Y(i-1))) ./ (1 + A);
end