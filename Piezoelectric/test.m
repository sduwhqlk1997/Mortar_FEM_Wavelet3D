%% 测试AssembleSAWMatFEM
clear
clc
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
ModelCoeff = 'D:\Code\M\Mortar_FEM_Wavelet\Piezoelectric\Data\ModelCoef2.mat';
load(ModelCoeff)
omega = 1e7;
c_LN = cell2mat(materials(2));
lamda_Al = cell2mat(materials(3));
mu_Al = cell2mat(materials(4));
e_LN = cell2mat(materials(5));
epcl_LN = cell2mat(materials(7));
density_LN = cell2mat(materials(8));
density_Al = cell2mat(materials(9));
ori1 = cell2mat(Geo(1));a1 = cell2mat(Geo(2));b1 = cell2mat(Geo(3));c1 = cell2mat(Geo(4));
ori2 = cell2mat(Geo(5));a2 = cell2mat(Geo(6));b2 = cell2mat(Geo(7));c2 = cell2mat(Geo(8));
Nx1 = 17; Ny1 = 2; Nz1 = 17;
Nx2 = 9; Nz2 = 5;

type="quadratic";
% type="linear";
[K,M,Dof_Index,~,~,~] = AssembleSAWMatFEM(c_LN,e_LN,epcl_LN,density_LN,...
    lamda_Al,mu_Al,density_Al,...
    ori1,a1,b1,c1,Nx1,Ny1,Nz1,...
    ori2,a2,b2,c2,Nx2,Nz2,true,type);
