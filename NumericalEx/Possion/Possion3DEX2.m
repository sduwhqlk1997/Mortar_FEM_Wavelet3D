%% Possion方程算例（做区域分解，并用mortar算法拼接区域）
clear
% clc
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet\FEM_FEM';
addpath(genpath(currentPath));
% 构造算例
interval=[0,1;0,1;0,2]; % 总区域
syms u(x,y,z)
% u(x,y,z) = sin(x+y+z);
u(x,y,z) = sin(x).*sin(y).*sin(z)...
    .*(x-interval(1,1)).*(x-interval(1,2))...
    .*(y-interval(2,1)).*(y-interval(2,2))...
    .*(z-interval(3,1)).*(z-interval(3,2));
f = -(diff(u,x,2)+diff(u,y,2)+diff(u,z,2));
u_exact = matlabFunction(u);
f = matlabFunction(f);
% 几何信息
% 下部区域
ori1=[0;0;0];
a1=1; b1=1; c1=1;
interval1 = [ori1,ori1+[a1;b1;c1]];
% 上部区域
ori2=[0;0;1];
a2=1; b2=1; c2=1;
interval2 = [ori2,ori2+[a2;b2;c2]];
% 网格信息
% 下部区域
Nx1=16; Ny1=16; Nz1=16;
xStep1=linspace(interval1(1,1),interval1(1,2),Nx1+1)';
yStep1=linspace(interval1(2,1),interval1(2,2),Ny1+1)';
zStep1=linspace(interval1(3,1),interval1(3,2),Nz1+1)';
% 上部区域
Nx2=6; Ny2=6; Nz2=4;
xStep2=linspace(interval2(1,1),interval2(1,2),Nx2+1)';
yStep2=linspace(interval2(2,1),interval2(2,2),Ny2+1)';
zStep2=linspace(interval2(3,1),interval2(3,2),Nz2+1)';
% 子区域网格剖分
% FEMBaseType='linear';
FEMBaseType='quadratic';
% 下部区域
[P1,T1,Pb1,Tb1] = genMesh3D(ori1,a1,b1,c1,[],[],[],FEMBaseType,xStep1,yStep1,zStep1);
[P2,T2,Pb2,Tb2] = genMesh3D(ori2,a2,b2,c2,[],[],[],FEMBaseType,xStep2,yStep2,zStep2);
viewMesh(P1,T1)
hold on
viewMesh(P2,T2)
% 参考单元信息计算
order=4;
[pt,w] = genRefGauss3DCube(order);
[phi,Dphi] = BaseFunOnRefGauss(order,FEMBaseType);
% 矩阵、载荷向量组装
% 下部区域
[Jacobi1,Jacobi_inv1,y1] = AffineToCuboid(pt,P1,T1,1);
A11 = AssembleFEMMatrix(1,Pb1,Tb1,[1,0,0],[1,0,0],phi,Dphi,Jacobi1,Jacobi_inv1,w);
A22 = AssembleFEMMatrix(1,Pb1,Tb1,[0,1,0],[0,1,0],phi,Dphi,Jacobi1,Jacobi_inv1,w);
A33 = AssembleFEMMatrix(1,Pb1,Tb1,[0,0,1],[0,0,1],phi,Dphi,Jacobi1,Jacobi_inv1,w);
K1 = A11+A22+A33;
F1 = AssembleFEMVector(f,Pb1,Tb1,phi,Jacobi1,y1,w);
% 上部区域
[Jacobi2,Jacobi_inv2,y2] = AffineToCuboid(pt,P2,T2,1);
A11 = AssembleFEMMatrix(1,Pb2,Tb2,[1,0,0],[1,0,0],phi,Dphi,Jacobi2,Jacobi_inv2,w);
A22 = AssembleFEMMatrix(1,Pb2,Tb2,[0,1,0],[0,1,0],phi,Dphi,Jacobi2,Jacobi_inv2,w);
A33 = AssembleFEMMatrix(1,Pb2,Tb2,[0,0,1],[0,0,1],phi,Dphi,Jacobi2,Jacobi_inv2,w);
K2 = A11+A22+A33;
F2 = AssembleFEMVector(f,Pb2,Tb2,phi,Jacobi2,y2,w);
K=blkdiag(K1,K2);
% 组装Lagrange乘子矩阵
MasterSiderInfo={interval1,FEMBaseType,xStep1,yStep1,zStep1,Pb1,Tb1};
SlaveSiderInfo={interval2,FEMBaseType,xStep2,yStep2,zStep2,Pb2,Tb2};
[MatMasterSider,MatSlaveSider] = ConnectFEMFEMbyMortar(MasterSiderInfo,SlaveSiderInfo,order);
B=[MatMasterSider;-MatSlaveSider];
% B=[-MatSlaveSider;MatMasterSider];
K=[K,B;...
    B.',sparse(size(B,2),size(B,2))];
F=[F1;F2;zeros(size(B,2),1)];
% 处理Dirichlet边界条件
Pb=[Pb1;Pb2];
boundary_index=Pb(:,1)==interval(1,1)|Pb(:,1)==interval(1,2)|...
    Pb(:,2)==interval(2,1)|Pb(:,2)==interval(2,2)|...
    Pb(:,3)==interval(3,1)|Pb(:,3)==interval(3,2);
F = F - sum((u_exact(Pb(boundary_index,1),Pb(boundary_index,2),...
    Pb(boundary_index,3))).'.*K(:,boundary_index),2);
K(boundary_index,:)=[];
K(:,boundary_index)=[];
F(boundary_index,:)=[];
Pb(boundary_index,:)=[];
boundary1=boundary_index(1:size(Pb1,1));
boundary2=boundary_index(size(Pb1,1)+1:end);
Pb1(boundary1,:)=[];
Pb2(boundary2,:)=[];
% 计算数值解
uh=K\F;
uh=uh(1:size(Pb,1));
err=norm(u_exact(Pb(:,1),Pb(:,2),Pb(:,3))-uh,inf)/norm(u_exact(Pb(:,1),Pb(:,2),Pb(:,3)),inf)