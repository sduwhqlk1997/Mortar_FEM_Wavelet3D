clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
% 区域设置
ori = [0;0;0];
a=2; b=1; c=1;
interval = [ori,ori+[a;b;c]];
syms u(x,y,z)
u(x,y,z) = sin(x+y+z).*(x-interval(1,1)).*(x-interval(1,2))...
    .*(y-interval(2,1)).*(y-interval(2,2))...
    .*(z-interval(3,1)).*(z-interval(3,2));
% u(x,y,z) = (sin(x)+sin(y)+sin(z)).*(x-interval(1,1)).*(x-interval(1,2))...
%     .*(y-interval(2,1)).*(y-interval(2,2))...
%     .*(z-interval(3,1)).*(z-interval(3,2));
% u(x,y,z) = (x.^2+y.^2+z.^2).*(x-interval(1,1)).*(x-interval(1,2))...
%     .*(y-interval(2,1)).*(y-interval(2,2))...
%     .*(z-interval(3,1)).*(z-interval(3,2));
f = -(diff(u,x,2)+diff(u,y,2)+diff(u,z,2));
u_exact = matlabFunction(u);
f = matlabFunction(f);
order=4;
WaveBaseType='quadratic';
j0=2;
% 下部区域
ori1=ori;
a1=a/2; b1=b; c1=c;
interval1=[ori1,ori1+[a1;b1;c1]];
J1=8;
jcoef1=[1;1;1]; % j1coef*j1+j2coef*j2+j3coef*j3<=J
[~,jcoef_min1]=min(jcoef1);
jcoef_other1=[1;2;3];
jcoef_other1(jcoef_min1)=[];
j1=J1-(jcoef1(jcoef_other1(1))+jcoef1(jcoef_other1(2)))*j0;
% j1=J1-2*j0;
% 上部区域
ori2=ori+[a1;0;0];
a2=a/2; b2=b; c2=c;
interval2=[ori2,ori2+[a2;b2;c2]];
J2=8;
jcoef2=[1;1;1]; % j1coef*j1+j2coef*j2+j3coef*j3<=J
[~,jcoef_min2]=min(jcoef2);
jcoef_other2=[1;1;1];
jcoef_other2(jcoef_min2)=[];
j2=J2-(jcoef2(jcoef_other2(1))+jcoef2(jcoef_other2(2)))*j0;
% 公共存储区
%  生成一维基函数的信息
tic
[BaseFunIndex_ref,IntegralBlock_ref] = GenWaveletBaseIndex1D([0,1],j0,...
    max([j1,j2]),WaveBaseType);
Table = GenWaveNonZeroIndex1D(BaseFunIndex_ref);
%  生成一维参考基本矩阵
K00_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,0,0,WaveBaseType,WaveBaseType,order);
K10_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,1,0,WaveBaseType,WaveBaseType,order);
K11_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,1,1,WaveBaseType,WaveBaseType,order);
toc
% 下部区域矩阵组装
[DofIndex_sparse1,Table_sparse1] = GenWaveletSparseInfo3D(j0,J1,...
    BaseFunIndex_ref,BaseFunIndex_ref,BaseFunIndex_ref,jcoef1);
% 生成三维参考基本矩阵
K00_ref=AssembleWaveMatrix3D(DofIndex_sparse1,Table_sparse1,K00_1D,K00_1D,K00_1D);
K11_ref=AssembleWaveMatrix3D(DofIndex_sparse1,Table_sparse1,K11_1D,K00_1D,K00_1D);
K12_ref=AssembleWaveMatrix3D(DofIndex_sparse1,Table_sparse1,K10_1D,K10_1D.',K00_1D);
K22_ref=AssembleWaveMatrix3D(DofIndex_sparse1,Table_sparse1,K00_1D,K11_1D,K00_1D);
K13_ref=AssembleWaveMatrix3D(DofIndex_sparse1,Table_sparse1,K10_1D,K00_1D,K10_1D.');
K23_ref=AssembleWaveMatrix3D(DofIndex_sparse1,Table_sparse1,K00_1D,K10_1D,K10_1D.');
K33_ref=AssembleWaveMatrix3D(DofIndex_sparse1,Table_sparse1,K00_1D,K00_1D,K11_1D);
%  区域剖分
Block = {[interval1(1,1),interval1(1,2)],...
    [interval1(2,1),interval1(2,2)],...
    [interval1(3,1),interval1(3,2)]};
%  生成总体基本矩阵
[~,K11,~,K22,~,~,K33,...
    Dof_index1,~,~,~] = ...
    CombineSomeWaveDomain3D(Block,...
    K00_ref,K11_ref,K12_ref,K22_ref,K13_ref,K23_ref,K33_ref,...
    DofIndex_sparse1,BaseFunIndex_ref);
K1 = K11+K22+K33;
toc
F1 = AssembleWaveVectorMatching3D(f,Dof_index1,WaveBaseType,order);
% 上部区域矩阵组装
[DofIndex_sparse2,Table_sparse2] = GenWaveletSparseInfo3D(j0,J2,...
    BaseFunIndex_ref,BaseFunIndex_ref,BaseFunIndex_ref,jcoef2);
% 生成三维参考基本矩阵
K00_ref=AssembleWaveMatrix3D(DofIndex_sparse2,Table_sparse2,K00_1D,K00_1D,K00_1D);
K11_ref=AssembleWaveMatrix3D(DofIndex_sparse2,Table_sparse2,K11_1D,K00_1D,K00_1D);
K12_ref=AssembleWaveMatrix3D(DofIndex_sparse2,Table_sparse2,K10_1D,K10_1D.',K00_1D);
K22_ref=AssembleWaveMatrix3D(DofIndex_sparse2,Table_sparse2,K00_1D,K11_1D,K00_1D);
K13_ref=AssembleWaveMatrix3D(DofIndex_sparse2,Table_sparse2,K10_1D,K00_1D,K10_1D.');
K23_ref=AssembleWaveMatrix3D(DofIndex_sparse2,Table_sparse2,K00_1D,K10_1D,K10_1D.');
K33_ref=AssembleWaveMatrix3D(DofIndex_sparse2,Table_sparse2,K00_1D,K00_1D,K11_1D);
%  区域剖分
Block = {[interval2(1,1),interval2(1,2)],...
    [interval2(2,1),interval2(2,2)],...
    [interval2(3,1),interval2(3,2)]};
%  生成总体基本矩阵
[~,K11,~,K22,~,~,K33,...
    Dof_index2,~,~,~] = ...
    CombineSomeWaveDomain3D(Block,...
    K00_ref,K11_ref,K12_ref,K22_ref,K13_ref,K23_ref,K33_ref,...
    DofIndex_sparse2,BaseFunIndex_ref);
K2 = K11+K22+K33;
F2 = AssembleWaveVectorMatching3D(f,Dof_index2,WaveBaseType,order);
% 处理两个区域的Dirichlet边界条件
%  下部区域
% boundary_index1 = ...
%     unique([ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,1);interval(2,1),interval(2,2);interval(3,:)],Dof_index1);...
%     ExtractWaveBoundaryIndex3D([interval(1,2),interval(1,2);interval(2,1),interval(2,2);interval(3,:)],Dof_index1);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,1);interval(3,:)],Dof_index1);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,2),interval(2,2);interval(3,:)],Dof_index1);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,2);interval(3,1),interval1(3,1)],Dof_index1)]);

% boundary_index1 = ...
%     unique([ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,1);interval(2,1),interval(2,2);interval(3,:)],Dof_index1);...
%     ExtractWaveBoundaryIndex3D([interval(1,2),interval(1,2);interval(2,1),interval(2,2);interval(3,:)],Dof_index1);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,1);interval(3,:)],Dof_index1);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,2);interval(3,1),interval(3,1)],Dof_index1);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,2);interval(3,2),interval1(3,2)],Dof_index1)]);

boundary_index1 = ...
    unique([ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,1);interval(2,1),interval(2,2);interval(3,:)],Dof_index1);...
    ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,1);interval(3,:)],Dof_index1);...
    ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,2),interval(2,2);interval(3,:)],Dof_index1);...
    ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,2);interval(3,1),interval(3,1)],Dof_index1);...
    ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,2);interval(3,2),interval1(3,2)],Dof_index1)]);

index_I1=(1:size(Dof_index1,1))';
index_I1(boundary_index1)=[];
K1(boundary_index1,:)=[];
K1(:,boundary_index1)=[];
F1(boundary_index1)=[];
Dof_index1(boundary_index1,:)=[];
%  上部区域
% boundary_index2 = ...
%     unique([ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,1);interval(2,1),interval(2,2);interval(3,:)],Dof_index2);...
%     ExtractWaveBoundaryIndex3D([interval(1,2),interval(1,2);interval(2,1),interval(2,2);interval(3,:)],Dof_index2);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,1);interval(3,:)],Dof_index2);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,2),interval(2,2);interval(3,:)],Dof_index2);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,2);interval(3,2),interval(3,2)],Dof_index2)]);
% boundary_index2 = ...
%     unique([ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,1);interval(2,1),interval(2,2);interval(3,:)],Dof_index2);...
%     ExtractWaveBoundaryIndex3D([interval(1,2),interval(1,2);interval(2,1),interval(2,2);interval(3,:)],Dof_index2);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,2),interval(2,2);interval(3,:)],Dof_index2);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,2);interval(3,1),interval(3,1)],Dof_index2);...
%     ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,2);interval(3,2),interval(3,2)],Dof_index2)]);

boundary_index2 = ...
    unique([ExtractWaveBoundaryIndex3D([interval(1,2),interval(1,2);interval(2,1),interval(2,2);interval(3,:)],Dof_index2);...
    ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,1);interval(3,:)],Dof_index2);...
    ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,2),interval(2,2);interval(3,:)],Dof_index2);...
    ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,2);interval(3,1),interval(3,1)],Dof_index2);...
    ExtractWaveBoundaryIndex3D([interval(1,1),interval(1,2);interval(2,1),interval(2,2);interval(3,2),interval(3,2)],Dof_index2)]);

index_I2=(1:size(Dof_index2,1))';
index_I2(boundary_index2)=[];
K2(boundary_index2,:)=[];
K2(:,boundary_index2)=[];
F2(boundary_index2)=[];
Dof_index2(boundary_index2,:)=[];
% 组装Lagrange乘子连接矩阵
MasterSideInfo={interval1,Dof_index1(:,4:6),Dof_index1(:,7:9),Dof_index1(:,1:3),Dof_index1(:,16:18),ones(size(Dof_index1,1),1)};
SlaveSideInfo={interval2,Dof_index2(:,4:6),Dof_index2(:,7:9),Dof_index2(:,1:3),Dof_index2(:,16:18),ones(size(Dof_index2,1),1)};
[MatMaster,MatSlave] = Connect2WaveByDualMortar(MasterSideInfo,SlaveSideInfo);
B=[MatMaster;-MatSlave];
% 组装总体矩阵
K=blkdiag(K1,K2);
K=[K,B;B',sparse(size(B,2),size(B,2))];
F=[F1;F2;zeros(size(B,2),1)];
% 计算数值解
uh=K\F;
uh1=uh(1:size(Dof_index1,1)); % 区域1
uh2=uh(size(Dof_index1,1)+1:size(Dof_index1,1)+size(Dof_index2,1)); % 区域2
% uh1=sortrows([[index_I1,uh1];[boundary_index1,zeros(length(boundary_index1),1)]]);
% uh1=uh1(:,2);
% uh2=sortrows([[index_I2,uh2];[boundary_index2,zeros(length(boundary_index2),1)]]);
% uh2=uh2(:,2);
% 计算误差
u1=@(x,y,z)ApproxWaveFun3D(x,y,z,uh1,Dof_index1,[0,0,0],WaveBaseType);
u2=@(x,y,z)ApproxWaveFun3D(x,y,z,uh2,Dof_index2,[0,0,0],WaveBaseType);
xx=interval1(1,1):(a1/50):interval1(1,2);
yy=interval1(2,1):(b1/50):interval1(2,2);
zz=interval1(3,1):(c1/50):interval1(3,2);
[xx,yy,zz] = meshgrid(xx,yy,zz);
xx=xx(:); yy=yy(:); zz=zz(:);
err1 = norm(u1(xx(:),yy(:),zz(:))-u_exact(xx(:),yy(:),zz(:)),inf)...
    /norm(u_exact(xx(:),yy(:),zz(:)),inf)

xx=interval2(1,1):(a2/50):interval2(1,2);
yy=interval2(2,1):(b2/50):interval2(2,2);
zz=interval2(3,1):(c2/50):interval2(3,2);
[xx,yy,zz] = meshgrid(xx,yy,zz);
xx=xx(:); yy=yy(:); zz=zz(:);
err2 = norm(u2(xx(:),yy(:),zz(:))-u_exact(xx(:),yy(:),zz(:)),inf)...
    /norm(u_exact(xx(:),yy(:),zz(:)),inf)