% 与FEM3D中的算例Possion3DRefSolForMortar.m
% 沿x方向划分三个区域
clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
% 中间区域上表面的边界条件
% u_M_top=@(x,y,z)2-((x-1.5).^2+(y-0.5).^2); % EX1
% u_M_top=@(x,y,z)sin(x+y); % EX2
u_M_top=@(x,y,z)ones(size(x,1),size(x,2)); % EX3
% 划分区域
%  左
ori_L=[0;0;0];
a_L=1;
b_L=1;
c_L=1;
interval_L=[ori_L,ori_L+[a_L;b_L;c_L]];
%  中
ori_M=[ori_L(1)+a_L;0;0];
a_M=1;
b_M=1;
c_M=1;
interval_M=[ori_M,ori_M+[a_M;b_M;c_M]];
%  右
ori_R=[ori_M(1)+a_M;0;0];
a_R=1;
b_R=1;
c_R=1;
interval_R=[ori_R,ori_R+[a_R;b_R;c_R]];
% 
order=4;
WaveBaseType='quadratic';
j0=2;
% 各区域的离散信息
J=7;
%  左
J_L=J;
jcoef_L=[1;1;1];
j_L = MaxJ2Sparse(jcoef_L,j0,J_L);
%  中
J_M=J;
jcoef_M=[1;1;1];
j_M = MaxJ2Sparse(jcoef_M,j0,J_M);
%  右
J_R=J;
jcoef_R=[1;1;1];
j_R = MaxJ2Sparse(jcoef_R,j0,J_R);
% 公共存储区
%  生成一维基函数信息
[BaseFunIndex_ref,IntegralBlock_ref] = GenWaveletBaseIndex1D([0,1],j0,...
    max([j_L,j_M,j_R]),WaveBaseType);
Table = GenWaveNonZeroIndex1D(BaseFunIndex_ref);
%  生成一维参考基本矩阵
K00_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,0,0,WaveBaseType,WaveBaseType,order);
K10_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,1,0,WaveBaseType,WaveBaseType,order);
K11_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,1,1,WaveBaseType,WaveBaseType,order);
% 组装各区域刚度矩阵
%  左
[DofIndex_sparse_L,Table_sparse_L] = GenWaveletSparseInfo3D(j0,J_L,...
    BaseFunIndex_ref,BaseFunIndex_ref,BaseFunIndex_ref,jcoef_L);
K00=AssembleWaveMatrix3D(DofIndex_sparse_L,Table_sparse_L,K00_1D,K00_1D,K00_1D);
K11=AssembleWaveMatrix3D(DofIndex_sparse_L,Table_sparse_L,K11_1D,K00_1D,K00_1D);
K12=AssembleWaveMatrix3D(DofIndex_sparse_L,Table_sparse_L,K10_1D,K10_1D.',K00_1D);
K22=AssembleWaveMatrix3D(DofIndex_sparse_L,Table_sparse_L,K00_1D,K11_1D,K00_1D);
K13=AssembleWaveMatrix3D(DofIndex_sparse_L,Table_sparse_L,K10_1D,K00_1D,K10_1D.');
K23=AssembleWaveMatrix3D(DofIndex_sparse_L,Table_sparse_L,K00_1D,K10_1D,K10_1D.');
K33=AssembleWaveMatrix3D(DofIndex_sparse_L,Table_sparse_L,K00_1D,K00_1D,K11_1D);
Block = {[interval_L(1,1),interval_L(1,2)],...
    [interval_L(2,1),interval_L(2,2)],...
    [interval_L(3,1),interval_L(3,2)]};
[~,K11,~,K22,~,~,K33,...
    Dof_index_L,~,~,~] = ...
    CombineSomeWaveDomain3D(Block,...
    K00,K11,K12,K22,K13,K23,K33,...
    DofIndex_sparse_L,BaseFunIndex_ref);
Dof_index_L(:,19:24)=[];
K_L=K11+K22+K33;
%  中
[DofIndex_sparse_M,Table_sparse_M] = GenWaveletSparseInfo3D(j0,J_M,...
    BaseFunIndex_ref,BaseFunIndex_ref,BaseFunIndex_ref,jcoef_M);
K00=AssembleWaveMatrix3D(DofIndex_sparse_M,Table_sparse_M,K00_1D,K00_1D,K00_1D);
K11=AssembleWaveMatrix3D(DofIndex_sparse_M,Table_sparse_M,K11_1D,K00_1D,K00_1D);
K12=AssembleWaveMatrix3D(DofIndex_sparse_M,Table_sparse_M,K10_1D,K10_1D.',K00_1D);
K22=AssembleWaveMatrix3D(DofIndex_sparse_M,Table_sparse_M,K00_1D,K11_1D,K00_1D);
K13=AssembleWaveMatrix3D(DofIndex_sparse_M,Table_sparse_M,K10_1D,K00_1D,K10_1D.');
K23=AssembleWaveMatrix3D(DofIndex_sparse_M,Table_sparse_M,K00_1D,K10_1D,K10_1D.');
K33=AssembleWaveMatrix3D(DofIndex_sparse_M,Table_sparse_M,K00_1D,K00_1D,K11_1D);
Block = {[interval_M(1,1),interval_M(1,2)],...
    [interval_M(2,1),interval_M(2,2)],...
    [interval_M(3,1),interval_M(3,2)]};
[~,K11,~,K22,~,~,K33,...
    Dof_index_M,~,~,~] = ...
    CombineSomeWaveDomain3D(Block,...
    K00,K11,K12,K22,K13,K23,K33,...
    DofIndex_sparse_M,BaseFunIndex_ref);
Dof_index_M(:,19:24)=[];
K_M=K11+K22+K33;
%  右
[DofIndex_sparse_R,Table_sparse_R] = GenWaveletSparseInfo3D(j0,J_R,...
    BaseFunIndex_ref,BaseFunIndex_ref,BaseFunIndex_ref,jcoef_R);
K00=AssembleWaveMatrix3D(DofIndex_sparse_R,Table_sparse_R,K00_1D,K00_1D,K00_1D);
K11=AssembleWaveMatrix3D(DofIndex_sparse_R,Table_sparse_R,K11_1D,K00_1D,K00_1D);
K12=AssembleWaveMatrix3D(DofIndex_sparse_R,Table_sparse_R,K10_1D,K10_1D.',K00_1D);
K22=AssembleWaveMatrix3D(DofIndex_sparse_R,Table_sparse_R,K00_1D,K11_1D,K00_1D);
K13=AssembleWaveMatrix3D(DofIndex_sparse_R,Table_sparse_R,K10_1D,K00_1D,K10_1D.');
K23=AssembleWaveMatrix3D(DofIndex_sparse_R,Table_sparse_R,K00_1D,K10_1D,K10_1D.');
K33=AssembleWaveMatrix3D(DofIndex_sparse_R,Table_sparse_R,K00_1D,K00_1D,K11_1D);
Block = {[interval_R(1,1),interval_R(1,2)],...
    [interval_R(2,1),interval_R(2,2)],...
    [interval_R(3,1),interval_R(3,2)]};
[~,K11,~,K22,~,~,K33,...
    Dof_index_R,~,~,~] = ...
    CombineSomeWaveDomain3D(Block,...
    K00,K11,K12,K22,K13,K23,K33,...
    DofIndex_sparse_R,BaseFunIndex_ref);
Dof_index_R(:,19:24)=[];
K_R=K11+K22+K33;

% 处理中部区域Diriclet边界条件
N_M=size(K_M,1);
boundary_index=find(Dof_index_M(:,18)==interval_M(3,2)&...
    Dof_index_M(:,16)>=interval_M(1,1) & Dof_index_M(:,16)<=interval_M(1,2));
[uB_M,indexI_M,indexB_M] = WaveTraceFun3D(u_M_top,boundary_index,Dof_index_M,WaveBaseType);
F_M=zeros(N_M,1);
F_M=F_M-sum(uB_M'.*K_M(:,boundary_index),2);
K_M(boundary_index,:)=[];
K_M(:,boundary_index)=[];
F_M(boundary_index)=[];
% 处理左边区域的右上角边
N_L=size(K_L,1);
[indexRT,indexML]=find(Dof_index_L(:,16)==Dof_index_M(indexB_M,16)'&...
    Dof_index_L(:,17)==Dof_index_M(indexB_M,17)'&...
    Dof_index_L(:,18)==Dof_index_M(indexB_M,18)');
uB_L=uB_M(indexML);
indexB_L=indexRT;
indexI_L=(1:N_L)';
indexI_L(indexB_L)=[];
F_L=zeros(N_L,1);
F_L=F_L-sum(uB_L'.*K_L(:,indexB_L),2);
K_L(indexB_L,:)=[];
K_L(:,indexB_L)=[];
F_L(indexB_L)=[];
% 处理右边区域的左上角边
N_R=size(K_R,1);
[indexLT,indexMR]=find(Dof_index_R(:,16)==Dof_index_M(indexB_M,16)'&...
    Dof_index_R(:,17)==Dof_index_M(indexB_M,17)'&...
    Dof_index_R(:,18)==Dof_index_M(indexB_M,18)');
uB_R=uB_M(indexMR);
indexB_R=indexLT;
indexI_R=(1:N_R)';
indexI_R(indexB_R)=[];
F_R=zeros(N_R,1);
F_R=F_R-sum(uB_R'.*K_R(:,indexB_R),2);
K_R(indexB_R,:)=[];
K_R(:,indexB_R)=[];
F_R(indexB_R)=[];
% 组装Lagrange乘子连接矩阵
LeftSideInfo={interval_L,Dof_index_L(indexI_L,4:6),Dof_index_L(indexI_L,7:9),...
    Dof_index_L(indexI_L,1:3),Dof_index_L(indexI_L,16:18),ones(length(indexI_L),1)};
MidSideInfo={interval_M,Dof_index_M(indexI_M,4:6),Dof_index_M(indexI_M,7:9),...
    Dof_index_M(indexI_M,1:3),Dof_index_M(indexI_M,16:18),ones(length(indexI_M),1)};
RightSideInfo={interval_R,Dof_index_R(indexI_R,4:6),Dof_index_R(indexI_R,7:9),...
    Dof_index_R(indexI_R,1:3),Dof_index_R(indexI_R,16:18),ones(length(indexI_R),1)};
[B_L,B_M1] = Connect2WaveByDualMortar(LeftSideInfo,MidSideInfo);
[B_R,B_M2] = Connect2WaveByDualMortar(RightSideInfo,MidSideInfo);
B=[B_L,zeros(size(B_L,1),size(B_M2,2));...
    -B_M1,B_M2;...
    zeros(size(B_R,1),size(B_M1,2)),-B_R];
K=blkdiag(K_L,K_M,K_R);
K=[K,B;...
    B',sparse(size(B,2),size(B,2))];
F=[F_L;F_M;F_R;zeros(size(B,2),1)];
% 计算数值解
uh=K\F;
uh_L=uh(1:length(indexI_L));
uh(1:length(indexI_L))=[];
uh_M=uh(1:length(indexI_M));
uh(1:length(indexI_M))=[];
uh_R=uh(1:length(indexI_R));
uh(1:length(indexI_R))=[];

uh_L=sortrows([indexI_L,uh_L;indexB_L,uB_L]);
uh_L=uh_L(:,2);

uh_R=sortrows([indexI_R,uh_R;indexB_R,uB_R]);
uh_R=uh_R(:,2);

uh_M=sortrows([indexI_M,uh_M;indexB_M,uB_M]);
uh_M=uh_M(:,2);

u_L=@(x,y,z)ApproxWaveFun3D(x,y,z,uh_L,Dof_index_L,[0,0,0],WaveBaseType);
u_M=@(x,y,z)ApproxWaveFun3D(x,y,z,uh_M,Dof_index_M,[0,0,0],WaveBaseType);
u_R=@(x,y,z)ApproxWaveFun3D(x,y,z,uh_R,Dof_index_R,[0,0,0],WaveBaseType);

% 绘制数值解图
figure(2)
DrawCubeDomain(interval_L)
hold on
DrawCubeDomain(interval_M)
DrawCubeDomain(interval_R)
load('D:\Code\M\Mortar_FEM_Wavelet\NumericalEx\Possion\Data\FEMRefSol3ReigionsEx3.mat')
% 左
index_FEM_L=Dof_index_FEM(:,2)>=interval_L(1,1)&Dof_index_FEM(:,2)<interval_L(1,2);
xx=Dof_index_FEM(index_FEM_L,2);
yy=Dof_index_FEM(index_FEM_L,3);
zz=Dof_index_FEM(index_FEM_L,4);
u_Wave=u_L(xx,yy,zz);
err_L=norm(u_Wave-uh_FEM(index_FEM_L),inf)/norm(uh_FEM(index_FEM_L),inf)
scatter3(xx,yy,zz,10,u_Wave,'filled')
% 中
index_FEM_M=Dof_index_FEM(:,2)>=interval_M(1,1)&Dof_index_FEM(:,2)<=interval_M(1,2);
xx=Dof_index_FEM(index_FEM_M,2);
yy=Dof_index_FEM(index_FEM_M,3);
zz=Dof_index_FEM(index_FEM_M,4);
u_Wave=u_M(xx,yy,zz);
err_M=norm(u_Wave-uh_FEM(index_FEM_M),inf)/norm(uh_FEM(index_FEM_M),inf)
scatter3(xx,yy,zz,10,u_Wave,'filled')
% 右
index_FEM_R=Dof_index_FEM(:,2)>interval_R(1,1)&Dof_index_FEM(:,2)<=interval_R(1,2);
xx=Dof_index_FEM(index_FEM_R,2);
yy=Dof_index_FEM(index_FEM_R,3);
zz=Dof_index_FEM(index_FEM_R,4);
u_Wave=u_R(xx,yy,zz);
err_R=norm(u_Wave-uh_FEM(index_FEM_R),inf)/norm(uh_FEM(index_FEM_R),inf)
scatter3(xx,yy,zz,10,u_Wave,'filled')
% 绘制数值解图
% figure(3)
% DrawCubeDomain(interval_L)
% hold on
% DrawCubeDomain(interval_M)
% DrawCubeDomain(interval_R)
% load('D:\Code\M\Mortar_FEM_Wavelet\NumericalEx\Possion\Data\FEMRefSol3ReigionsEx2.mat')
% % 左
% index_FEM_L=Dof_index_FEM(:,2)>=interval_L(1,1)&Dof_index_FEM(:,2)<interval_L(1,2);
% xx=Dof_index_FEM(index_FEM_L,2);
% yy=Dof_index_FEM(index_FEM_L,3);
% zz=Dof_index_FEM(index_FEM_L,4);
% u_Wave=u_L(xx,yy,zz);
% u_Wave=u_Wave-uh_FEM(index_FEM_L);
% scatter3(xx,yy,zz,10,u_Wave,'filled')
% % 中
% index_FEM_M=Dof_index_FEM(:,2)>=interval_M(1,1)&Dof_index_FEM(:,2)<=interval_M(1,2);
% xx=Dof_index_FEM(index_FEM_M,2);
% yy=Dof_index_FEM(index_FEM_M,3);
% zz=Dof_index_FEM(index_FEM_M,4);
% u_Wave=u_M(xx,yy,zz);
% u_Wave=u_Wave-uh_FEM(index_FEM_M);
% scatter3(xx,yy,zz,10,u_Wave,'filled')
% % 右
% index_FEM_R=Dof_index_FEM(:,2)>interval_R(1,1)&Dof_index_FEM(:,2)<=interval_R(1,2);
% xx=Dof_index_FEM(index_FEM_R,2);
% yy=Dof_index_FEM(index_FEM_R,3);
% zz=Dof_index_FEM(index_FEM_R,4);
% u_Wave=u_R(xx,yy,zz);
% u_Wave=u_Wave-uh_FEM(index_FEM_R);
% scatter3(xx,yy,zz,10,u_Wave,'filled')


