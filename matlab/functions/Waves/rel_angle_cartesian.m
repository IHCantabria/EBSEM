function relD=rel_angle_cartesian(waveD, batiD)

%     ###########################################################################    
%     # Relative angle (in degrees) between wave direction and bathymetry with 
%     # angles in cartesian coordinates, angle between [180,-180], 
%     # 0 is in EAST and positive counterclockwise.
%     #
%     # INPUT:
%     # waveD:    wave angle in Cartesian notation.
%     # batiD:    bathymetry angle (normal to the shoreline) in Cartesian notation.
%     #
%     # OUTPUT:
%     # relD:     relative wave angle between wave and bathymetry, 0 is the bathymetry and positive counterclockwise.
%     ###########################################################################    
% caso = 0;
% if isinstance(waveD,float) or isinstance(waveD,int):
%     waveD = np.asarray([waveD])
%     caso = 1
% elseif isinstance(waveD,list)
%     waveD = np.asarray(waveD)
%     caso = 2;
% end

relD = waveD - batiD;
relD(relD>180)=relD(relD>180)-360;
relD(relD<-180)=relD(relD<-180)+360;

% if caso == 1
%     relD = relD.tolist()[0]
% elseif caso == 2
%     relD = relD.tolist()
% end
   
%     return relD
end