clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
ModelCoeff = 'D:\Code\M\Mortar_FEM_Wavelet\Piezoelectric\Data\ModelCoef2.mat';
load(ModelCoeff,'Geo')
ori1 = cell2mat(Geo(1));a1 = cell2mat(Geo(2));b1 = cell2mat(Geo(3));c1 = cell2mat(Geo(4));
ori2 = cell2mat(Geo(5));a2 = cell2mat(Geo(6));b2 = cell2mat(Geo(7));c2 = cell2mat(Geo(8));
% b1=100*1e-6;
% b2=b1;
interval_sub = [ori1.',ori1.'+[a1;b1;c1]];
interval_ele = [ori2.',ori2.'+[a2;b2;c2]];
DrawCubeDomain(interval_sub)
DrawCubeDomain(interval_ele)
% [P_sub,T_sub] = genMesh3D(ori1,a1,b1,c1,2,1001,2,"quadratic");
% [P_ele,T_ele] = genMesh3D(ori2,a2,b2,c2,2,1001,2,"quadratic");
% viewMesh(P_sub,T_sub)
% viewMesh(P_ele,T_ele)