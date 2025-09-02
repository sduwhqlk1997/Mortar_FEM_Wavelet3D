function [K,M,Dof_Index] = AssembleEleMatFEM(lamda_Al,mu_Al,density_Al,...
    ori2,a2,b2,c2,Nx2,Ny2,Nz2,...
    xStep_ele,yStep_ele,zStep_ele,type)
%ASSEMBLEELEMATFEM 组装电极（等势体）区域的刚度矩阵，只有线弹性部分
%   type：基函数类型
%   Dof_Index：第一列自由度类型1:u,2:v,3:w,4:phi；2-4列坐标
%% 公共存储区
interval=[ori2(1),ori2(1)+a2;ori2(2),ori2(2)+b2;ori2(3),ori2(3)+c2];
if ~exist('xStep_ele','var')||isempty(xStep_ele)
    xStep_ele=linspace(interval(1,1),interval(1,2),Nx1);
end
if ~exist('yStep_ele','var')||isempty(yStep_ele)
    yStep_ele=linspace(interval(2,1),interval(2,2),Ny1);
end
if ~exist('zStep_ele','var')||isempty(zStep_ele)
    zStep_ele=linspace(interval(3,1),interval(3,2),Nz1);
end
[P_ele,T_ele,Pb_ele,Tb_ele] =...
    genMesh3D(ori2,a2,b2,c2,Nx2,Ny2,Nz2,type,xStep_ele,yStep_ele,zStep_ele);
switch type
    case "linear"
        load('D:\Code\M\Mortar_FEM_Wavelet\FEM\Data\LinearRefEleGaussOrder4.mat');
    case "quadratic"
        load('D:\Code\M\Mortar_FEM_Wavelet\FEM\Data\QuadraticRefEleGaussOrder4.mat');
end
[Jacobi_ele,Jacobi_inv_ele] = AffineToCuboid(pt,P_ele,T_ele,0);
%% 组装刚度矩阵
B00 = AssembleFEMMatrix(1,Pb_ele,Tb_ele,[0,0,0],[0,0,0],phi,Dphi,Jacobi_ele,Jacobi_inv_ele,w);
B11 = AssembleFEMMatrix(1,Pb_ele,Tb_ele,[1,0,0],[1,0,0],phi,Dphi,Jacobi_ele,Jacobi_inv_ele,w);
B12 = AssembleFEMMatrix(1,Pb_ele,Tb_ele,[1,0,0],[0,1,0],phi,Dphi,Jacobi_ele,Jacobi_inv_ele,w);
B13 = AssembleFEMMatrix(1,Pb_ele,Tb_ele,[1,0,0],[0,0,1],phi,Dphi,Jacobi_ele,Jacobi_inv_ele,w);
B22 = AssembleFEMMatrix(1,Pb_ele,Tb_ele,[0,1,0],[0,1,0],phi,Dphi,Jacobi_ele,Jacobi_inv_ele,w);
B23 = AssembleFEMMatrix(1,Pb_ele,Tb_ele,[0,1,0],[0,0,1],phi,Dphi,Jacobi_ele,Jacobi_inv_ele,w);
B33 = AssembleFEMMatrix(1,Pb_ele,Tb_ele,[0,0,1],[0,0,1],phi,Dphi,Jacobi_ele,Jacobi_inv_ele,w);
%   组装Muu
Muu_e = density_Al*B00;
M = blkdiag(Muu_e,Muu_e,Muu_e);
%   组装Kuu
mu11 = (lamda_Al+2*mu_Al)*B11+mu_Al*B22+mu_Al*B33;
mu12 = lamda_Al*B12.'+mu_Al*B12;
mu13 = lamda_Al*B13.'+mu_Al*B13;
mu22 = (lamda_Al+2*mu_Al)*B22+mu_Al*B11+mu_Al*B33;
mu23 = lamda_Al*B23.'+mu_Al*B23;
mu33 = (lamda_Al+2*mu_Al)*B33+mu_Al*B22+mu_Al*B11;
K = [mu11,mu12,mu13;mu12.',mu22,mu23;mu13.',mu23.',mu33];
%   生成自由度索引
Dof_Index = [kron([1;2;3],ones(size(Pb_ele,1),1)),kron([1;1;1],Pb_ele)];
end

