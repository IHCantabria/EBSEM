function [H2, DIR2, h2] = BreakingPropagation(H1,T1,DIR1,h1,ANGbati)
% BreakingPropagation(H1,T1,DIR1,h1,ANGbati,breakType)
%     ###########################################################################    
%     # Propagation of waves using linear theory assuming rectilinear and parallel bathymetry
%     #    
%     # INPUT:
%     # H1:        wave height.
%     # T1:        wave period.
%     # DIR1:      wave direction. Nautical convention.
%     # h1:        depth of wave conditions.
%     # ANGbati:   bathymetry angle, the normal of the shoreline. Cartesian notation
%     # breakType: type of breaking condition. Spectral or monochromatic.
%     #    
%     # OUTPUT:
%     # H2:        wave height during breaking. Wave period is assumed invariant due to linear theory
%     # DIR2:      wave direction during breaking. Nautical convention.
%     # h2:        depth of breaking
%     ###########################################################################    

% breaksT = {'spectral':0.45, 'mono':0.78};

% Bcoef=breaksT[breakType.lower()]
    
Bcoef=.55;
DIRrel = rel_angle_cartesian(nauticalDir2cartesianDir(DIR1),ANGbati);


h2l0 = H1./Bcoef; % # initial condition for breaking depth
    

H2 = zeros(numel(H1),1); 
DIR2 = zeros(numel(DIR1),1); 
h2 = zeros(numel(H1),1);


H2(h2l0>=h1) = H1(h2l0>=h1); 
DIR2(h2l0>=h1) = DIR1(h2l0>=h1); 
h2(h2l0>=h1) = h2l0(h2l0>=h1);   % # check that the initial depth is deeper than the breaking value

H2(H1<=0.1) = H1(H1<=0.1); 
DIR2(H1<=0.1) = DIR1(H1<=0.1);
h2(H1<=0.1) = h2l0(H1<=0.1);   % # check that the initial depth is deeper than the breaking value

propProf = find((abs(DIRrel)<=90) & (H1>0.1) & (h2l0<h1));
% #    print 'init'    
% #    print len(propProf)
% #    print 'end'
if ~isempty(propProf)
    myFun= @(x) LinearShoalBreak_Residual(x, H1(propProf), ...
        T1(propProf), DIR1(propProf), h1, ANGbati, ...
        Bcoef);
    h2l = JFNK(myFun,h2l0(propProf), 1e-3, 500);
    [~, H2l, DIR2l] = LinearShoalBreak_Residual(h2l, H1(propProf),...
        T1(propProf), DIR1(propProf), h1, ANGbati, Bcoef);                
    H2(propProf) = H2l;
    DIR2(propProf) = DIR2l;
    h2(propProf) = h2l;
end
    
%     return H2, DIR2, h2

end