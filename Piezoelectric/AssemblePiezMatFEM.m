function [K,M,F,Dof_Index] =...
    AssemblePiezMatFEM(c_LN,e_LN,epcl_LN,density_LN,kappa_bar_sq,...
    ori1,a1,b1,c1,Nx1,Ny1,Nz1,...
    xStep_sub,yStep_sub,zStep_sub,type,...
    RightHand)
%ASSEMBLEPIEZMATFEM 用FEM基组装压电基底区域的刚度矩阵
%   kappa_bar_sq:机电耦合系数，无量纲化方程或许需要，不做无量纲化忽略即可
%   Nx1,Ny1,Nz1：三个坐标方向的网格点数
%   type：基函数类型
%   Dof_Index：第一列自由度类型1:u,2:v,3:w,4:phi；2-4列坐标
%   RightHand：右端载荷，为4x1Cell，依次为u,v,w,phi的载荷
interval=[ori1(1),ori1(1)+a1;ori1(2),ori1(2)+b1;ori1(3),ori1(3)+c1];
if ~exist('xStep_sub','var')||isempty(xStep_sub)
    xStep_sub=linspace(interval(1,1),interval(1,2),Nx1);
end
if ~exist('yStep_sub','var')||isempty(yStep_sub)
    yStep_sub=linspace(interval(2,1),interval(2,2),Ny1);
end
if ~exist('zStep_sub','var')||isempty(zStep_sub)
    zStep_sub=linspace(interval(3,1),interval(3,2),Nz1);
end
[P_sub,T_sub,Pb_sub,Tb_sub] =...
    genMesh3D(ori1,a1,b1,c1,Nx1,Ny1,Nz1,type,xStep_sub,yStep_sub,zStep_sub);
%% 公共存储区
switch type
    case "linear"
        % load('D:\Code\M\Mortar_FEM_Wavelet\FEM\Data\LinearRefEleGaussOrder4.mat');
        load('LinearRefEleGaussOrder4.mat');
    case "quadratic"
        % load('D:\Code\M\Mortar_FEM_Wavelet\FEM\Data\QuadraticRefEleGaussOrder4.mat');
        load('QuadraticRefEleGaussOrder4.mat')
end
[Jacobi_sub,Jacobi_inv_sub,y] = AffineToCuboid(pt,P_sub,T_sub,1);
%% 组装刚度矩阵
A00 = AssembleFEMMatrix(1,Pb_sub,Tb_sub,[0,0,0],[0,0,0],phi,Dphi,Jacobi_sub,Jacobi_inv_sub,w);
A11 = AssembleFEMMatrix(1,Pb_sub,Tb_sub,[1,0,0],[1,0,0],phi,Dphi,Jacobi_sub,Jacobi_inv_sub,w);
A12 = AssembleFEMMatrix(1,Pb_sub,Tb_sub,[1,0,0],[0,1,0],phi,Dphi,Jacobi_sub,Jacobi_inv_sub,w);
A13 = AssembleFEMMatrix(1,Pb_sub,Tb_sub,[1,0,0],[0,0,1],phi,Dphi,Jacobi_sub,Jacobi_inv_sub,w);
A22 = AssembleFEMMatrix(1,Pb_sub,Tb_sub,[0,1,0],[0,1,0],phi,Dphi,Jacobi_sub,Jacobi_inv_sub,w);
A23 = AssembleFEMMatrix(1,Pb_sub,Tb_sub,[0,1,0],[0,0,1],phi,Dphi,Jacobi_sub,Jacobi_inv_sub,w);
A33 = AssembleFEMMatrix(1,Pb_sub,Tb_sub,[0,0,1],[0,0,1],phi,Dphi,Jacobi_sub,Jacobi_inv_sub,w);

%   组装Muu
Muu = density_LN*A00;
Muu = blkdiag(Muu,Muu,Muu);
%   组装Kuu
u11 = C(c_LN,1,1,1,1)*A11+C(c_LN,1,1,1,2)*(A12+A12.')+C(c_LN,1,1,1,3)*(A13+A13.')+ ...
    C(c_LN,1,2,1,2)*A22+C(c_LN,1,2,1,3)*(A23+A23.')+ ...
    C(c_LN,1,3,1,3)*A33;
u12 = C(c_LN,1,1,1,2)*A11+C(c_LN,1,1,2,2)*A12.'+C(c_LN,1,1,2,3)*A13.'+ ...
    C(c_LN,1,2,1,2)*A12+C(c_LN,1,2,2,2)*A22+C(c_LN,1,2,2,3)*A23.'+ ...
    C(c_LN,1,3,1,2)*A13+C(c_LN,1,3,2,2)*A23+C(c_LN,1,3,2,3)*A33;
u13 = C(c_LN,1,1,1,3)*A11+C(c_LN,1,1,2,3)*A12.'+C(c_LN,1,1,3,3)*A13.'+ ...
    C(c_LN,1,2,1,3)*A12+C(c_LN,1,2,2,3)*A22+C(c_LN,1,2,3,3)*A23.'+ ...
    C(c_LN,1,3,1,3)*A13+C(c_LN,1,3,2,3)*A23+C(c_LN,1,3,3,3)*A33;
u22 = C(c_LN,1,2,1,2)*A11+C(c_LN,1,2,2,2)*(A12+A12.')+C(c_LN,1,2,2,3)*(A13+A13.')+ ...
    C(c_LN,2,2,2,2)*A22+C(c_LN,2,2,2,3)*(A23+A23.')+ ...
    C(c_LN,2,3,2,3)*A33;
u23 = C(c_LN,1,2,1,3)*A11+C(c_LN,1,2,2,3)*A12.'+C(c_LN,1,2,3,3)*A13.'+ ...
    C(c_LN,2,2,1,3)*A12+C(c_LN,2,2,2,3)*A22+C(c_LN,2,2,3,3)*A23.'+ ...
    C(c_LN,2,3,1,3)*A13+C(c_LN,2,3,2,3)*A23+C(c_LN,2,3,3,3)*A33;
u33 = C(c_LN,1,3,1,3)*A11+C(c_LN,1,3,2,3)*(A12+A12.')+C(c_LN,1,3,3,3)*(A13+A13.')+ ...
    C(c_LN,2,3,2,3)*A22+C(c_LN,2,3,3,3)*(A23+A23.')+ ...
    C(c_LN,3,3,3,3)*A33;
Kuu = [u11,u12,u13;u12.',u22,u23;u13.',u23.',u33];
%   组装Kup
phi11 = E(e_LN,1,1,1)*A11+E(e_LN,2,1,1)*A12.'+E(e_LN,3,1,1)*A13.'+...
    E(e_LN,1,2,1)*A12+E(e_LN,2,2,1)*A22+E(e_LN,3,2,1)*A23.'+...
    E(e_LN,1,3,1)*A13+E(e_LN,2,3,1)*A23+E(e_LN,3,3,1)*A33;
phi21 = E(e_LN,1,1,2)*A11+E(e_LN,2,1,2)*A12.'+E(e_LN,3,1,2)*A13.'+...
    E(e_LN,1,2,2)*A12+E(e_LN,2,2,2)*A22+E(e_LN,3,2,2)*A23.'+...
    E(e_LN,1,3,2)*A13+E(e_LN,2,3,2)*A23+E(e_LN,3,3,2)*A33;
phi31 = E(e_LN,1,1,3)*A11+E(e_LN,2,1,3)*A12.'+E(e_LN,3,1,3)*A13.'+...
    E(e_LN,1,2,3)*A12+E(e_LN,2,2,3)*A22+E(e_LN,3,2,3)*A23.'+...
    E(e_LN,1,3,3)*A13+E(e_LN,2,3,3)*A23+E(e_LN,3,3,3)*A33;
Kup = [phi11;phi21;phi31];
%   组装Kp
Kp = epcl_LN(1,1)*A11+epcl_LN(1,2)*A12+epcl_LN(1,3)*A13+...
    epcl_LN(2,1)*A12.'+epcl_LN(2,2)*A22+epcl_LN(2,3)*A23+...
    epcl_LN(3,1)*A13.'+epcl_LN(3,2)*A23.'+epcl_LN(3,3)*A33;
%   形成总体刚度矩阵
if ~exist('kappa_bar_sq','var')||isempty(kappa_bar_sq)
    kappa_bar_sq=1;
end
K = [Kuu,kappa_bar_sq*Kup;...
    Kup.',-Kp];
%   形成总体质量阵
M = blkdiag(Muu,sparse(size(Kp,1),size(Kp,1)));
%   形成自由度索引
Dof_Index = [kron([1;2;3;4],ones(size(Pb_sub,1),1)),kron([1;1;1;1],Pb_sub)];
if exist('RightHand','var')&&~isempty(RightHand)
    f1=cell2mat(RightHand(1));
    f2=cell2mat(RightHand(2));
    f3=cell2mat(RightHand(3));
    f4=cell2mat(RightHand(4));

    F1 = AssembleFEMVector(f1,Pb_sub,Tb_sub,phi,Jacobi_sub,y,w);
    F2 = AssembleFEMVector(f2,Pb_sub,Tb_sub,phi,Jacobi_sub,y,w);
    F3 = AssembleFEMVector(f3,Pb_sub,Tb_sub,phi,Jacobi_sub,y,w);
    F4 = AssembleFEMVector(f4,Pb_sub,Tb_sub,phi,Jacobi_sub,y,w);

    F=[F1;F2;F3;F4];
else
    F=[];
end
end

