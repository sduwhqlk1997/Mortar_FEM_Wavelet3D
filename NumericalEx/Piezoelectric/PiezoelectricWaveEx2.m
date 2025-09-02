% 纯压电算例，不做区域分解
clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
% 上表面局部加1V电压
load('ModelCoef.mat')
% 压电基底几何参数
ori_p=cell2mat(Geo(1))';
a1=cell2mat(Geo(2));
b1=cell2mat(Geo(3));
c1=cell2mat(Geo(4));
% 上表面施加电压
V0=1;
% 区域
interval_p=[ori_p,ori_p+[a1;b1;c1]];
% 材料参数
omega_p=cell2mat(materials(1));
% omega_p=0;
c_p=cell2mat(materials(2));
e_p=cell2mat(materials(5));
epcl_p=cell2mat(materials(7)).*cell2mat(materials(6));
rho_p=cell2mat(materials(8));
MaterialPara={c_p,e_p,epcl_p,rho_p};
% 小波离散信息
order=4;
WaveBaseType='quadratic';
j0=2;
J=8;
jcoef=[1;1;1];
j=MaxJ2Sparse(jcoef,j0,J);

% 公共信息
[BaseFunIndex_ref,IntegralBlock_ref] = GenWaveletBaseIndex1D([0,1],j0,...
    j,WaveBaseType);
Table = GenWaveNonZeroIndex1D(BaseFunIndex_ref);

Mat1DRef={AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,0,0,WaveBaseType,WaveBaseType,order),...
    AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,1,0,WaveBaseType,WaveBaseType,order),...
    AssembleWaveMatrix1D(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
    IntegralBlock_ref,1,1,WaveBaseType,WaveBaseType,order)};
% Mat1DRef={AssembleWaveMatrix1D_2(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
%     0,0,WaveBaseType,WaveBaseType,order),...
%     AssembleWaveMatrix1D_2(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
%     1,0,WaveBaseType,WaveBaseType,order),...
%     AssembleWaveMatrix1D_2(@(x) 1,[0,1],Table,BaseFunIndex_ref,...
%     1,1,WaveBaseType,WaveBaseType,order)};
% 组装刚度矩阵
[K,M,BaseIndex,DofIndex] = AssemblePiezMatWave(Mat1DRef,BaseFunIndex_ref,...
    jcoef,J,j0,MaterialPara,interval_p);
% 处理Dirichlet边界条件
K=K-omega_p^2*M;
BottomIndex=find(BaseIndex(:,18)==interval_p(3,1));
TopIndex=find(BaseIndex(:,18)==interval_p(3,2));
DiriInfo={1,BottomIndex,0,0;2,BottomIndex,0,0;3,BottomIndex,0,0;4,BottomIndex,0,0;...
    4,TopIndex,@(x,y,z)ones(size(x,1),size(x,2))*V0,1};
[K,F,uIndexI,vIndexI,wIndexI,phiIndexI,uB,vB,wB,phiB] ...
    = TreatPiezDiriBoundWave3D(K,zeros(size(K,1),1),DiriInfo,BaseIndex,WaveBaseType);
% 计算数值解
uh=K\F;
% 生成各变量的数值解
u1=sortrows([uB;uIndexI(:,2),uh(uIndexI(:,1))]); u1=u1(:,2);
u2=sortrows([vB;vIndexI(:,2),uh(vIndexI(:,1))]); u2=u2(:,2);
u3=sortrows([wB;wIndexI(:,2),uh(wIndexI(:,1))]); u3=u3(:,2);
phih=sortrows([phiB;phiIndexI(:,2),uh(phiIndexI(:,1))]); phih=phih(:,2);
u1=@(x,y,z) ApproxWaveFun3D(x,y,z,u1,BaseIndex,[0,0,0],WaveBaseType);
u2=@(x,y,z) ApproxWaveFun3D(x,y,z,u2,BaseIndex,[0,0,0],WaveBaseType);
u3=@(x,y,z) ApproxWaveFun3D(x,y,z,u3,BaseIndex,[0,0,0],WaveBaseType);
phih=@(x,y,z) ApproxWaveFun3D(x,y,z,phih,BaseIndex,[0,0,0],WaveBaseType);
% 绘制数值解图
%  生成绘制点
% [~,~,Pb,~] = genMesh3D(ori_p,a1,b1,c1,...
%     cell2mat(Geo(9)),5,cell2mat(Geo(11)),"quadratic");
load('FEMSOL.mat')
xx=Table_FEM(:,1); yy=Table_FEM(:,2); zz=Table_FEM(:,3);
figure(1)
DrawCubeDomain(interval_p);
% Myfindnode(Pb)
hold on
u1=u1(xx,yy,zz); v1=u2(xx,yy,zz); w1=u3(xx,yy,zz);
Pos=[xx,yy,zz]+[u1,v1,w1]*5e3;
Disp=sqrt(u1.^2+v1.^2+w1.^2);
scatter3(Pos(:,1),Pos(:,2),Pos(:,3),10,Disp,'filled')
err1=norm(Table_FEM(:,4)-u1,inf)/norm(Table_FEM(:,4),inf)
err2=norm(Table_FEM(:,5)-v1,inf)/norm(Table_FEM(:,5),inf)
err3=norm(Table_FEM(:,6)-w1,inf)/norm(Table_FEM(:,6),inf)
xx=Table_FEM_phi(:,1); yy=Table_FEM_phi(:,2); zz=Table_FEM_phi(:,3);
phih=phih(xx,yy,zz);
errphi=norm(Table_FEM_phi(:,4)-phih,inf)/norm(Table_FEM_phi(:,4),inf)
% err1=norm(Table_FEM(:,4)-u1,inf)
% err2=norm(Table_FEM(:,5)-v1,inf)
% err3=norm(Table_FEM(:,6)-w1,inf)
% errall=norm(Table_FEM(:,4:6)-[u1,v1,w1],inf)
%.........未完待续


