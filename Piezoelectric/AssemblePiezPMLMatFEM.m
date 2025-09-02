function [K,M,Dof_index] = AssemblePiezPMLMatFEM(Df1,Df2,Df3,...
    c_pml,e_pml,epcl_pml,density_pml,...
    ori,a,b,c,Nx,Ny,Nz,...
    xStep_pml,yStep_pml,zStep_pml,type)
%ASSEMBLEPIEZPMLMATFEM 组装压电方程的PML区域矩阵
%   Df1,Df2,Df3:三个方向的PML函数的导数(注：Df_i为x_i的单变量函数)，依次为df_1(x_1)/dx_1、df_2(x_2)/dx_2和df_3(x_3)/dx_3
%   type：基函数类型
%   Dof_Index：第一列自由度类型1:u,2:v,3:w,4:phi；2-4列坐标
%% 生成网格步长
interval=[ori(1),ori(1)+a;ori(2),ori(2)+b;ori(3),ori(3)+c];
if ~exist('xStep_pml','var')||isempty(xStep_pml)
    xStep_pml=linspace(interval(1,1),interval(1,2),Nx);
end
if ~exist('yStep_pml','var')||isempty(yStep_pml)
    yStep_pml=linspace(interval(2,1),interval(2,2),Ny);
end
if ~exist('zStep_pml','var')||isempty(zStep_pml)
    zStep_pml=linspace(interval(3,1),interval(3,2),Nz);
end
%% 计算坐标变换的Jacobi矩阵的对角元
alpha1 = @(x,y,z) 1-1j*Df1(x);
alpha2 = @(x,y,z) 1-1j*Df2(y);
alpha3 = @(x,y,z) 1-1j*Df3(z);
%% 网格剖分
[P_pml,T_pml,Pb_pml,Tb_pml] = ...
    genMesh3D(ori,a,b,c,Nx,Ny,Nz,type,xStep_pml,yStep_pml,zStep_pml);
%% 公共存储区
switch type
    case "linear"
        load('D:\Code\M\Mortar_FEM_Wavelet\FEM\Data\LinearRefEleGaussOrder4.mat');
    case "quadratic"
        load('D:\Code\M\Mortar_FEM_Wavelet\FEM\Data\QuadraticRefEleGaussOrder4.mat');
end
[Jacobi_pml,Jacobi_inv_pml,y] = AffineToCuboid(pt,P_pml,T_pml,1);
%   计算所需非线性系数
coef00 = @(x,y,z) alpha1(x,y,z).*alpha2(x,y,z).*alpha3(x,y,z);
coef11 = @(x,y,z) alpha2(x,y,z).*alpha3(x,y,z)./alpha1(x,y,z);
coef12 = @(x,y,z) alpha3(x,y,z);
coef13 = @(x,y,z) alpha2(x,y,z);
coef22 = @(x,y,z) alpha1(x,y,z).*alpha3(x,y,z)./alpha2(x,y,z);
coef23 = @(x,y,z) alpha1(x,y,z);
coef33 = @(x,y,z) alpha1(x,y,z).*alpha2(x,y,z)./alpha3(x,y,z);
%   计算所需非线性系数在所有Gauss求积节点的值
Coef00 = reshape(coef00(reshape(y(:,:,1),[],1),reshape(y(:,:,2),[],1),reshape(y(:,:,3),[],1)),size(Tb_pml,1),[]);
Coef11 = reshape(coef11(reshape(y(:,:,1),[],1),reshape(y(:,:,2),[],1),reshape(y(:,:,3),[],1)),size(Tb_pml,1),[]);
Coef12 = reshape(coef12(reshape(y(:,:,1),[],1),reshape(y(:,:,2),[],1),reshape(y(:,:,3),[],1)),size(Tb_pml,1),[]);
Coef13 = reshape(coef13(reshape(y(:,:,1),[],1),reshape(y(:,:,2),[],1),reshape(y(:,:,3),[],1)),size(Tb_pml,1),[]);
Coef22 = reshape(coef22(reshape(y(:,:,1),[],1),reshape(y(:,:,2),[],1),reshape(y(:,:,3),[],1)),size(Tb_pml,1),[]);
Coef23 = reshape(coef23(reshape(y(:,:,1),[],1),reshape(y(:,:,2),[],1),reshape(y(:,:,3),[],1)),size(Tb_pml,1),[]);
Coef33 = reshape(coef33(reshape(y(:,:,1),[],1),reshape(y(:,:,2),[],1),reshape(y(:,:,3),[],1)),size(Tb_pml,1),[]);
%% 组装基本矩阵
A00 = AssembleFEMMatrix(Coef00,Pb_pml,Tb_pml,[0,0,0],[0,0,0],phi,Dphi,Jacobi_pml,Jacobi_inv_pml,w);
A11 = AssembleFEMMatrix(Coef11,Pb_pml,Tb_pml,[1,0,0],[1,0,0],phi,Dphi,Jacobi_pml,Jacobi_inv_pml,w);
A12 = AssembleFEMMatrix(Coef12,Pb_pml,Tb_pml,[1,0,0],[0,1,0],phi,Dphi,Jacobi_pml,Jacobi_inv_pml,w);
A13 = AssembleFEMMatrix(Coef13,Pb_pml,Tb_pml,[1,0,0],[0,0,1],phi,Dphi,Jacobi_pml,Jacobi_inv_pml,w);
A22 = AssembleFEMMatrix(Coef22,Pb_pml,Tb_pml,[0,1,0],[0,1,0],phi,Dphi,Jacobi_pml,Jacobi_inv_pml,w);
A23 = AssembleFEMMatrix(Coef23,Pb_pml,Tb_pml,[0,1,0],[0,0,1],phi,Dphi,Jacobi_pml,Jacobi_inv_pml,w);
A33 = AssembleFEMMatrix(Coef33,Pb_pml,Tb_pml,[0,0,1],[0,0,1],phi,Dphi,Jacobi_pml,Jacobi_inv_pml,w);
%   组装Muu
Muu=density_pml*A00;
Muu=blkdiag(Muu,Muu,Muu);
%   组装Kuu
u11 = C(c_pml,1,1,1,1)*A11+C(c_pml,1,1,1,2)*(A12+A12.')+C(c_pml,1,1,1,3)*(A13+A13.')+ ...
    C(c_pml,1,2,1,2)*A22+C(c_pml,1,2,1,3)*(A23+A23.')+ ...
    C(c_pml,1,3,1,3)*A33;
u12 = C(c_pml,1,1,1,2)*A11+C(c_pml,1,1,2,2)*A12.'+C(c_pml,1,1,2,3)*A13.'+ ...
    C(c_pml,1,2,1,2)*A12+C(c_pml,1,2,2,2)*A22+C(c_pml,1,2,2,3)*A23.'+ ...
    C(c_pml,1,3,1,2)*A13+C(c_pml,1,3,2,2)*A23+C(c_pml,1,3,2,3)*A33;
u13 = C(c_pml,1,1,1,3)*A11+C(c_pml,1,1,2,3)*A12.'+C(c_pml,1,1,3,3)*A13.'+ ...
    C(c_pml,1,2,1,3)*A12+C(c_pml,1,2,2,3)*A22+C(c_pml,1,2,3,3)*A23.'+ ...
    C(c_pml,1,3,1,3)*A13+C(c_pml,1,3,2,3)*A23+C(c_pml,1,3,3,3)*A33;
u22 = C(c_pml,1,2,1,2)*A11+C(c_pml,1,2,2,2)*(A12+A12.')+C(c_pml,1,2,2,3)*(A13+A13.')+ ...
    C(c_pml,2,2,2,2)*A22+C(c_pml,2,2,2,3)*(A23+A23.')+ ...
    C(c_pml,2,3,2,3)*A33;
u23 = C(c_pml,1,2,1,3)*A11+C(c_pml,1,2,2,3)*A12.'+C(c_pml,1,2,3,3)*A13.'+ ...
    C(c_pml,2,2,1,3)*A12+C(c_pml,2,2,2,3)*A22+C(c_pml,2,2,3,3)*A23.'+ ...
    C(c_pml,2,3,1,3)*A13+C(c_pml,2,3,2,3)*A23+C(c_pml,2,3,3,3)*A33;
u33 = C(c_pml,1,3,1,3)*A11+C(c_pml,1,3,2,3)*(A12+A12.')+C(c_pml,1,3,3,3)*(A13+A13.')+ ...
    C(c_pml,2,3,2,3)*A22+C(c_pml,2,3,3,3)*(A23+A23.')+ ...
    C(c_pml,3,3,3,3)*A33;
Kuu = [u11,u12,u13;u12.',u22,u23;u13.',u23.',u33];
%   组装Kup
phi11 = E(e_pml,1,1,1)*A11+E(e_pml,2,1,1)*A12.'+E(e_pml,3,1,1)*A13.'+...
    E(e_pml,1,2,1)*A12+E(e_pml,2,2,1)*A22+E(e_pml,3,2,1)*A23.'+...
    E(e_pml,1,3,1)*A13+E(e_pml,2,3,1)*A23+E(e_pml,3,3,1)*A33;
phi21 = E(e_pml,1,1,2)*A11+E(e_pml,2,1,2)*A12.'+E(e_pml,3,1,2)*A13.'+...
    E(e_pml,1,2,2)*A12+E(e_pml,2,2,2)*A22+E(e_pml,3,2,2)*A23.'+...
    E(e_pml,1,3,2)*A13+E(e_pml,2,3,2)*A23+E(e_pml,3,3,2)*A33;
phi31 = E(e_pml,1,1,3)*A11+E(e_pml,2,1,3)*A12.'+E(e_pml,3,1,3)*A13.'+...
    E(e_pml,1,2,3)*A12+E(e_pml,2,2,3)*A22+E(e_pml,3,2,3)*A23.'+...
    E(e_pml,1,3,3)*A13+E(e_pml,2,3,3)*A23+E(e_pml,3,3,3)*A33;
Kup = [phi11;phi21;phi31];
%   组装Kp
Kp = epcl_pml(1,1)*A11+epcl_pml(1,2)*A12+epcl_pml(1,3)*A13+...
    epcl_pml(2,1)*A12.'+epcl_pml(2,2)*A22+epcl_pml(2,3)*A23+...
    epcl_pml(3,1)*A13.'+epcl_pml(3,2)*A23.'+epcl_pml(3,3)*A33;
%% 组装总体刚度矩阵和质量矩阵
K = [Kuu,Kup;...
    Kup.',-Kp];
M = blkdiag(Muu,0.*Kp);
%% 自由度索引
Dof_index = [kron([1;2;3;4],ones(size(Pb_pml,1),1)),kron([1;1;1;1],Pb_pml)];
end

