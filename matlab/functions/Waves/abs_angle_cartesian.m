function waveD = abs_angle_cartesian(relD, batiD)
%     ###########################################################################    
%     # Absolute angle in cartesian notation, angle between [180,-180], 
%     # 0 is in EAST and positive counterclockwise.
%     # From a relative angle from wave and bathymetry.
%     # The same as rel_angle_cartesian(relD,-1*batiD)
%     # INPUT:
%     # relD:     relative wave angle between wave and bathymetry, 0 is the bathymetry and positive counterclockwise.
%     # batiD:    bathymetry angle (normal to the shoreline) in Cartesian notation.
%     #
%     # OUTPUT:
%     # waveD:    wave angle in Cartesian notation.
%     ###########################################################################    

waveD = relD + batiD;
waveD(waveD>180)=waveD(waveD>180)-360;
waveD(waveD<-180)=waveD(waveD<-180)+360;
    
end