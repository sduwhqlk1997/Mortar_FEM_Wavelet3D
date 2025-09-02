function [xStep_sub,yStep_sub,zStep_sub,xStep_ele,yStep_ele,zStep_ele] = genMeshStepOfSAW(ori1,a1,b1,c1,Nx1,Ny1,Nz1,ori2,a2,c2,Nx2,Nz2)
%GENMESHSTEPOFSAW 生成SAW器件的网格步长
%   输入SAW器件的几何尺寸，输出基底和电极的网格步长向量
% y方向
yStep_sub=linspace(ori1(2),ori1(2)+b1,Ny1);
yStep_ele = yStep_sub;
% z方向
zStep_sub = linspace(ori1(3),ori1(3)+c1,Nz1);
zStep_ele = linspace(ori2(3),ori2(3)+c2,Nz2);
% x方向
x1 = linspace(ori1(1),ori2(1),floor((Nx1-Nx2)/2)+1);
x2 = linspace(ori2(1),ori2(1)+a2,Nx2);
xStep_ele = x2;
x3 = linspace(ori2(1)+a2,ori1(1)+a1,ceil((Nx1-Nx2)/2)+1);
xStep_sub = [x1,x2(2:end-1),x3];
end

