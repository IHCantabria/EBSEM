function nDir = cartesianDir2nauticalDir(cDir)
%     ###########################################################################    
%     # Cartesian convention with 0 in East and positive counterclockwise TO
%     # Nautical convention with 0 in North and positive clockwise. 
%     ###########################################################################    

nDir = 90.-cDir;
nDir(nDir<0) = 360+nDir(nDir<0);

end