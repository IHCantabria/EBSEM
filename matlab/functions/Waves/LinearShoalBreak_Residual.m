function [res, H2l, DIR2l]=LinearShoalBreak_Residual(h2l, H1, T1, DIR1, h1, ANGbati, Bcoef)

[H2l,DIR2l]= LinearShoal(H1, T1, DIR1, h1, h2l, ANGbati);
H2comp = h2l.*Bcoef;
res = H2l-H2comp;
% res=res';
end