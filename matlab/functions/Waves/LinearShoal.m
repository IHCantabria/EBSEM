function [H2, DIR2]=LinearShoal(H1, T1, DIR1, h1, h2, ANGbati)
%     ###########################################################################    
%     # Wave shoaling and refraction applying linear theory with parallel, rectilinear bathymetry.
%     #    
%     # INPUT:
%     # H1:        initial wave height.
%     # T1:        wave period.
%     # DIR1:      initial wave direction. Nautical convention.
%     # h1:        initial depth of wave conditions.
%     # h2:        final depth of wave conditions.
%     # ANGbati:   bathymetry angle, the normal of the shoreline. Cartesian convention
%     #
%     # OUTPUT:
%     # H2:        wave height during breaking. Wave period is assumed invariant due to linear theory.
%     # DIR2:      wave direction during breaking. Nautical convention.
%     ###########################################################################
relDir1 = rel_angle_cartesian(nauticalDir2cartesianDir(DIR1),ANGbati);
L1=hunt(h1,T1);
L2=hunt(h2,T1);
CG1 = GroupCelerity(L1,T1,h1);
CG2 = GroupCelerity(L2,T1,h2);
relDir2 = Snell_Law(L1,L2,relDir1);
KS = sqrt(CG1./CG2);
KR = sqrt(cos(relDir1.*pi./180.)./cos(relDir2.*pi./180.));
H2 = H1.*KS.*KR;
DIR2 = cartesianDir2nauticalDir(abs_angle_cartesian(relDir2,ANGbati));

end