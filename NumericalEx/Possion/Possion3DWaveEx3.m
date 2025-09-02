% 不做区域分解，所有面具有Dirichlet边界
clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
% 区域设置
ori = [0;0;0];
a=1; b=1; c=1;
interval=[ori(1)-a,ori(1)+a;ori(2)-b,ori(2)+b;ori(3)-c,ori(3)+c];

syms x y z
u(x,y,z) = sin(x+y+z);
g(x,y,z) = -(diff(u,x,2)+diff(u,y,2)+diff(u,z,2));
u_exact = matlabFunction(u);
f = matlabFunction(g);
% 设置一维基函数
J=7;
j0=2; j=J-2*j0;
BaseType="quadratic";
% 生成一维基函数信息
[BaseFunIndex_ref,IntegralBlock_ref] = GenWaveletBaseIndex1D([0,1],j0,j,BaseType);
Table = GenWaveNonZeroIndex1D(BaseFunIndex_ref);
% 生成一维参考基本矩阵
order=4;
K00_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,0,0,BaseType,BaseType,order);
K10_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,1,0,BaseType,BaseType,order);
K11_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,1,1,BaseType,BaseType,order);
% 生成三维稀疏基信息
[DofIndex_sparse,Table_sparse] = GenWaveletSparseInfo3D(j0,J,...
    BaseFunIndex_ref,BaseFunIndex_ref,BaseFunIndex_ref,[1,1,1]);
% 生成三维参考基本矩阵
K00=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K00_1D,K00_1D,K00_1D);
K11=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K11_1D,K00_1D,K00_1D);
K12=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K10_1D,K10_1D.',K00_1D);
K22=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K00_1D,K11_1D,K00_1D);
K13=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K10_1D,K00_1D,K10_1D.');
K23=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K00_1D,K10_1D,K10_1D.');
K33=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K00_1D,K00_1D,K11_1D);
% 将基本矩阵变换到要求的区间
Block = {[interval(1,1),interval(1,2)],...
    [interval(2,1),interval(2,2)],...
    [interval(3,1),interval(3,2)]};
[K00,K11,K12,K22,K13,K23,K33,...
    Dof_index,P,T,Block_interval] = ...
    CombineSomeWaveDomain3D(Block,...
    K00,K11,K12,K22,K13,K23,K33,...
    DofIndex_sparse,BaseFunIndex_ref);
% 组装总体刚度矩阵
K = K11+K22+K33;
F = AssembleWaveVectorMatching3D(f,Dof_index,BaseType,order);
% 处理边界条件
boundary_index=find(Dof_index(:,16)==interval(1,1)|...
    Dof_index(:,16)==interval(1,2)|...
    Dof_index(:,17)==interval(2,1)|...
    Dof_index(:,17)==interval(2,2)|...
    Dof_index(:,18)==interval(3,1)|...
    Dof_index(:,18)==interval(3,2));
Dof_index(:,19:24)=[];
[u_B,index_I,index_B] = WaveTraceFun3D(u_exact,boundary_index,Dof_index,BaseType);
F = F-sum(u_B'.*K(:,boundary_index),2);
K(boundary_index,:)=[];
K(:,boundary_index)=[];
F(boundary_index)=[];
% 求解线性方程组
u_I = K\F;
u_num=sortrows([[index_I,u_I];[boundary_index,u_B]]);
u_num = u_num(:,2);
% 生成数值解
u_h = @(x,y,z)ApproxWaveFun3D(x,y,z,u_num,Dof_index,[0,0,0],BaseType);
xx=interval(1,1):(a/20):interval(1,2);
yy=interval(2,1):(b/20):interval(2,2);
zz=interval(3,1):(c/20):interval(3,2);
[xx,yy,zz] = meshgrid(xx,yy,zz);
err = norm(u_h(xx(:),yy(:),zz(:))-u_exact(xx(:),yy(:),zz(:)),inf)...
    /norm(u_exact(xx(:),yy(:),zz(:)),inf)

