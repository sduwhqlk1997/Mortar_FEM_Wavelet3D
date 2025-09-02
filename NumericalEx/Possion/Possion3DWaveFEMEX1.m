% 上部区域用FEM离散，下部区域用SWGM离散，使用Mortar方法连接两个区域
clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
% 构造算例
ori1 = [0;0;0];
a1 = 1;b1=1;c1=2;
interval = [ori1,ori1+[a1;b1;c1]];
syms u(x,y,z)
% u(x,y,z) = sin(x+y+z);
u(x,y,z) = (sin(x)+sin(y)+sin(z)).*(x-interval(1,1)).*(x-interval(1,2))...
    .*(y-interval(2,1)).*(y-interval(2,2))...
    .*(z-interval(3,1)).*(z-interval(3,2));
% u(x,y,z) = (x.^2+y.^2+z.^2).*(x-interval(1,1)).*(x-interval(1,2))...
%     .*(y-interval(2,1)).*(y-interval(2,2))...
%     .*(z-interval(3,1)).*(z-interval(3,2));
f = -(diff(u,x,2)+diff(u,y,2)+diff(u,z,2));
u_exact = matlabFunction(u);
f = matlabFunction(f);
order=4;
% 有限元部分区域信息
FEMBaseType='quadratic';
oriFEM=ori1+[0;0;c1/2];
aFEM=a1; bFEM=b1; cFEM=c1/2;
intervalFEM=[oriFEM,oriFEM+[aFEM;bFEM;cFEM]];
Nx = 16; Ny=16; Nz=8;
% 小波部分区域信息
WaveBaseType='quadratic';
oriWave=ori1;
aWave=a1; bWave=b1; cWave=c1/2;
intervalWave=[oriWave,oriWave+[aWave;bWave;cWave]];
J=8;
j0=2;
j=J-2*j0;
% j=J;
% 有限元区域矩阵组装
xStepFEM=linspace(intervalFEM(1,1),intervalFEM(1,2),Nx+1)';
yStepFEM=linspace(intervalFEM(2,1),intervalFEM(2,2),Ny+1)';
zStepFEM=linspace(intervalFEM(3,1),intervalFEM(3,2),Nz+1)';
[P,T,Pb,Tb] = genMesh3D(oriFEM,aFEM,bFEM,cFEM,[],[],[],FEMBaseType,xStepFEM,yStepFEM,zStepFEM);
[pt,w] = genRefGauss3DCube(order);
[phi,Dphi] = BaseFunOnRefGauss(order,FEMBaseType);
[Jacobi,Jacobi_inv,y] = AffineToCuboid(pt,P,T,1);
A11 = AssembleFEMMatrix(1,Pb,Tb,[1,0,0],[1,0,0],phi,Dphi,Jacobi,Jacobi_inv,w);
A22 = AssembleFEMMatrix(1,Pb,Tb,[0,1,0],[0,1,0],phi,Dphi,Jacobi,Jacobi_inv,w);
A33 = AssembleFEMMatrix(1,Pb,Tb,[0,0,1],[0,0,1],phi,Dphi,Jacobi,Jacobi_inv,w);
KFEM = A11+A22+A33;
FFEM = AssembleFEMVector(f,Pb,Tb,phi,Jacobi,y,w);
% Dof_indexFEM = [ones(size(Pb,1),1),Pb,(1:size(Pb,1)).'];
% 小波区域矩阵组装
[BaseFunIndexWave_ref,~] = GenWaveletBaseIndex1D([0,1],j0,j,WaveBaseType);
TableWave = GenWaveNonZeroIndex1D(BaseFunIndexWave_ref); % 一维基函数信息
K00_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],TableWave,BaseFunIndexWave_ref,0,0,...一维参考矩阵
    WaveBaseType,WaveBaseType,order);
K10_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],TableWave,BaseFunIndexWave_ref,1,0,...
    WaveBaseType,WaveBaseType,order);
K11_1D = AssembleWaveMatrix1D(@(x) 1,[0,1],TableWave,BaseFunIndexWave_ref,1,1,...
    WaveBaseType,WaveBaseType,order);
[DofIndexWave_sparse,TableWave_sparse] = GenWaveletSparseInfo3D(j0,J,...稀疏基信息
    BaseFunIndexWave_ref,BaseFunIndexWave_ref,BaseFunIndexWave_ref);
K00_ref=AssembleWaveMatrix3D(DofIndexWave_sparse,TableWave_sparse,K00_1D,K00_1D,K00_1D); % 三维参考矩阵
K11_ref=AssembleWaveMatrix3D(DofIndexWave_sparse,TableWave_sparse,K11_1D,K00_1D,K00_1D);
K12_ref=AssembleWaveMatrix3D(DofIndexWave_sparse,TableWave_sparse,K10_1D,K10_1D.',K00_1D);
K22_ref=AssembleWaveMatrix3D(DofIndexWave_sparse,TableWave_sparse,K00_1D,K11_1D,K00_1D);
K13_ref=AssembleWaveMatrix3D(DofIndexWave_sparse,TableWave_sparse,K10_1D,K00_1D,K10_1D.');
K23_ref=AssembleWaveMatrix3D(DofIndexWave_sparse,TableWave_sparse,K00_1D,K10_1D,K10_1D.');
K33_ref=AssembleWaveMatrix3D(DofIndexWave_sparse,TableWave_sparse,K00_1D,K00_1D,K11_1D);
BlockWave = {[interval(1,1),interval(1,2)],... 小波元区域剖分
    [interval(2,1),interval(2,2)],...
    [interval(3,1),interval(3,2)]};
[~,K11,~,K22,~,~,K33,...
    Dof_indexWave,~,~,~] = ... 生成总体刚度矩阵
    CombineSomeWaveDomain3D(BlockWave,...
    K00_ref,K11_ref,K12_ref,K22_ref,K13_ref,K23_ref,K33_ref,...
    DofIndexWave_sparse,BaseFunIndexWave_ref);
KWave = K11+K22+K33;
FWave = AssembleWaveVectorMatching3D(f,Dof_indexWave,WaveBaseType,order);
% 生成Lagrange乘子矩阵
[BaseFunIndexx,IntegralBlockx]=GenWaveletBaseIndex1D(intervalWave(1,:),j0,j,WaveBaseType);
[BaseFunIndexy,IntegralBlocky]=GenWaveletBaseIndex1D(intervalWave(2,:),j0,j,WaveBaseType);
[BaseFunIndexz,IntegralBlockz]=GenWaveletBaseIndex1D(intervalWave(3,:),j0,j,WaveBaseType);
WaveSideInfo={intervalWave,WaveBaseType,BaseFunIndexx,BaseFunIndexy,BaseFunIndexz,...
    DofIndexWave_sparse,IntegralBlockx,IntegralBlocky,IntegralBlockz};
FEMSideInfo={intervalFEM,FEMBaseType,xStepFEM,yStepFEM,zStepFEM,Pb,Tb};
[MatWaveSide,MatFEMSide] = ConnectWaveletFEMbyMortar(WaveSideInfo,FEMSideInfo,order);
% 组装总体矩阵
K=blkdiag(KWave,KFEM);
B=[MatWaveSide;-MatFEMSide];
K=[K,B;B.',sparse(size(B,2),size(B,2))];
F=[FWave;FFEM;zeros(size(B,2),1)];
% 处理齐次Dirichlet边界条件
Pb_total=[Dof_indexWave(:,16:18),zeros(size(Dof_indexWave,1),1);Pb,ones(size(Pb,1),1)];
boundary_index=Pb_total(:,1)==interval(1,1)|Pb_total(:,1)==interval(1,2)|...
    Pb_total(:,2)==interval(2,1)|Pb_total(:,2)==interval(2,2)|...
    Pb_total(:,3)==interval(3,1)|Pb_total(:,3)==interval(3,2);
boundaryWavelet=boundary_index(1:size(Dof_indexWave,1));
boundaryFEM=boundary_index(size(Dof_indexWave,1)+1:end);

K(boundary_index,:)=[];
K(:,boundary_index)=[];
F(boundary_index)=[];
Pb_total(boundary_index,:)=[];
Dof_indexWave(boundaryWavelet,:)=[];
Pb(boundaryFEM,:)=[];

% 计算数值解
uh=K\F;
uh=uh(1:size(Pb_total,1));
uhWave=uh(Pb_total(:,end)==0);
uhFEM=uh(Pb_total(:,end)==1);
uWave = @(x,y,z)ApproxWaveFun3D(x,y,z,uhWave,Dof_indexWave,[0,0,0],WaveBaseType);% 小波区域
% 提取有限元交界面处的格点
% 交界面处误差
indexInterface=Pb(:,3)==1;
PbInterface=Pb(indexInterface,:);
uhFEMInterface=uhFEM(indexInterface,:);
uhWAVEInterface=uWave(PbInterface(:,1),PbInterface(:,2),PbInterface(:,3));
ExactInterface=u_exact(PbInterface(:,1),PbInterface(:,2),PbInterface(:,3));
errInterface=norm(uhFEMInterface-uhWAVEInterface,inf)/norm(uhFEMInterface,inf);
% 小波区域的误差
xx=(intervalWave(1,1):0.1:intervalWave(1,2))';
yy=(intervalWave(2,1):0.1:intervalWave(2,2))';
zz=(intervalWave(3,1):0.1:intervalWave(3,2))';
uNumWave=uWave(xx,yy,zz);
ExactWave=u_exact(xx,yy,zz);
errWave=norm(uNumWave-ExactWave,inf)/norm(ExactWave,inf)
% 有限元区域的误差
ExactFEM=u_exact(Pb(:,1),Pb(:,2),Pb(:,3));
errFEM=norm(uhFEM-ExactFEM,inf)/norm(ExactFEM,inf)
% PbFEM=Pb(Pb(:,end)==1,:);
% err=norm(uhFEM-u_exact(PbFEM(:,1),PbFEM(:,2),PbFEM(:,3)),inf)/...
%     norm(u_exact(PbFEM(:,1),PbFEM(:,2),PbFEM(:,3)),inf); % 有限元部分误差

