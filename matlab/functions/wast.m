function wsf = wast(hb,D50)

%    ###########################################################################    
%    # Width of the active surf zone
%    # hb:   depth of closure
%    # D50:  mean sediment grain size (m)
%    ###########################################################################    
%   #    hb = hb+CM ????? see why is introducing the tidal range in the width of the active surf zone    

wsf=(hb./ADEAN(D50)).^(3./2.);

end