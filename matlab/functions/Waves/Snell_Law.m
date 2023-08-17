function alpha = Snell_Law(L1,L2,alpha1)
%     ###########################################################################    
%     # Wave refraction using snell law.
%     #    
%     # INPUT:
%     # L1:     initial wave length.
%     # L1:     final wave length.
%     # alpha1: initial wave dir. Cartesian notation.
%     #
%     # OUTPUT:
%     # alpha1: final wave dir. Cartesian notation.
%     ###########################################################################    
alpha=asin(L2.*sin(alpha1.*pi./180.)./L1).*180./pi;

end