function Cg = GroupCelerity(L,T,h)
%     ###########################################################################    
%     # CELERITY GROUP
%     # L: wave lenght.
%     # T: wave period.
%     # h: depth of wave conditions.
%     ###########################################################################       

c = L./T;
k = 2 .*pi./L;
N = 1+2 .*k.*h./sinh(2 .*k.*h);
Cg = c /2 .* N;

end