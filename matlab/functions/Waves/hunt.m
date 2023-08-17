function [L,L0]=hunt(d,T)
%CALCULO DE LA LONGITUD DDE ONDA POR EL METODO DE HUNT

%% constantes

g=9.81; %[m/s^2]

%% Calculos

L0=g.*T.^2./(2*pi);

G=2.*pi.*(d./L0);

p=[.067,.0864,.4622,.6522,1];

F= G+1./polyval(p,G);

L=T.*(g.*d./F).^.5;

end