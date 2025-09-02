%% Possion方程算例（不做区域分解）
clear
% clc
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
% 构造算例
ori1 = [0,0,0];
a1 = 1;b1=1;c1=2;
interval = [ori1.',ori1.'+[a1;b1;c1]];
syms u(x,y,z)
% u(x,y,z) = sin(x+y+z);
u(x,y,z) = 10*sin(x+y+z).*(x-interval(1,1)).*(x-interval(1,2))...
    .*(y-interval(2,1)).*(y-interval(2,2))...
    .*(z-interval(3,1)).*(z-interval(3,2));
f = -(diff(u,x,2)+diff(u,y,2)+diff(u,z,2));
u_exact = matlabFunction(u);
f = matlabFunction(f);

% 几何信息
% ori1 = [0,0,0];
% a1 = 1;b1=1;c1=1;
order=4;
FEMBaseType='linear';
% FEMBaseType="quadratic";
N = 16;
xStep=linspace(interval(1,1),interval(1,2),N+1)';
yStep=linspace(interval(2,1),interval(2,2),N+1)';
zStep=linspace(interval(3,1),interval(3,2),N+1)';
[P,T,Pb,Tb] = genMesh3D(ori1,a1,b1,c1,[],[],[],FEMBaseType,xStep,yStep,zStep);
% Nx1=N+1;Ny1=N+1;Nz1=N+1;
% [P,T,Pb,Tb] = genMesh3D(ori1,a1,b1,c1,Nx1,Ny1,Nz1,FEMBaseType);

[pt,w] = genRefGauss3DCube(order);
[phi,Dphi] = BaseFunOnRefGauss(order,FEMBaseType);
% load('QuadraticRefEleGaussOrder4.mat')

[Jacobi,Jacobi_inv,y] = AffineToCuboid(pt,P,T,1);
% 组装矩阵
A11 = AssembleFEMMatrix(1,Pb,Tb,[1,0,0],[1,0,0],phi,Dphi,Jacobi,Jacobi_inv,w);
A22 = AssembleFEMMatrix(1,Pb,Tb,[0,1,0],[0,1,0],phi,Dphi,Jacobi,Jacobi_inv,w);
A33 = AssembleFEMMatrix(1,Pb,Tb,[0,0,1],[0,0,1],phi,Dphi,Jacobi,Jacobi_inv,w);
K = A11+A22+A33;
F = AssembleFEMVector(f,Pb,Tb,phi,Jacobi,y,w);
Dof_index = [ones(size(Pb,1),1),Pb,(1:size(Pb,1)).'];
% 处理边界条件
boundary_index = find(Dof_index(:,2)==ori1(1)|Dof_index(:,2)==ori1(1)+a1|...
    Dof_index(:,3)==ori1(2)|Dof_index(:,3)==ori1(2)+b1|...
    Dof_index(:,4)==ori1(3)|Dof_index(:,4)==ori1(3)+c1);
% F = F - sum((u_exact(Dof_index(boundary_index,2),Dof_index(boundary_index,3),...
%     Dof_index(boundary_index,4))).'.*K(:,boundary_index),2);
K(boundary_index,:) = [];
K(:,boundary_index) = [];
F(boundary_index) = [];
Bound_Info = {ones(size(boundary_index,1),1),Dof_index(boundary_index,2:4),...
    u_exact(Dof_index(boundary_index,2),Dof_index(boundary_index,3),...
    Dof_index(boundary_index,4)),boundary_index};
Dof_index(boundary_index,:)=[];
% 求解方程
u_h = K\F;
[Dof_index,u_h] = CompleteNumerResult(Dof_index,u_h,Bound_Info);
% 生成数值解
err = ...
    norm(u_exact(Dof_index(:,2),Dof_index(:,3),Dof_index(:,4))-u_h,inf)/norm(u_exact(Dof_index(:,2),Dof_index(:,3),Dof_index(:,4)),inf)