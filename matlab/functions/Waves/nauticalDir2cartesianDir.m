function cDir=nauticalDir2cartesianDir(nDir)

%     ###########################################################################    
%     # Nautical convention with 0 in North and positive clockwise TO 
%     # Cartesian convention with 0 in East and positive counterclockwise.
%     ###########################################################################    
% caso = 0;
% if isinstance(nDir,float) || isinstance(nDir,int)
%     nDir = np.asarray([nDir])
%     caso = 1
% elseif isinstance(nDir,list)
%     nDir = np.asarray(nDir)
%     caso = 2
% end
        
cDir = 90-nDir;
cDir(cDir<(-180)) = 360+cDir(cDir<(-180));
% if caso == 1
%     cDir = cDir.tolist()[0]
% elseif caso == 2
%     cDir = cDir.tolist()
% end
        
%     return cDir


end