% 有界域压电方程算例
clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
% 上表面局部加1V电压
load('ModelCoef.mat')
% 压电基底几何参数
ori_p=cell2mat(Geo(1))';
a_p=cell2mat(Geo(2));
b_p=cell2mat(Geo(3));
c_p=cell2mat(Geo(4));
% 施加电压的区域
V0=1;
ori_e=cell2mat(Geo(5))';
a_e=cell2mat(Geo(6));
b_e=cell2mat(Geo(7));
% 划分区域
%  左
interval_L=[ori_p(1),ori_e(1);...
    ori_p(2),ori_p(2)+b_p;...
    ori_p(3),ori_p(3)+c_p];
%  中（电）
interval_E=[ori_e(1),ori_e(1)+a_e;...
    interval_L(2:3,:)];
%  右
interval_R=[interval_E(1,2),ori_p(1)+a_p;...
    interval_L(2:3,:)];
% 材料参数
omega_p=cell2mat(materials(1));
c_p=cell2mat(materials(2));
e_p=cell2mat(materials(5));
epcl_p=cell2mat(materials(7)).*cell2mat(materials(6));
rho_p=cell2mat(materials(8));
MaterialPara={c_p,e_p,epcl_p,rho_p};
% 小波离散信息
order=4;
WaveBaseType='quadratic';
j0=2;
J_L=7; J_E=7; J_R=7;
jLcoef=[1;1;1]; jEcoef=[1;1;1]; jRcoef=[1;1;1];
j_L=MaxJ2Sparse(jLcoef,j0,J_L);
j_E=MaxJ2Sparse(jEcoef,j0,J_E);
j_R=MaxJ2Sparse(jRcoef,j0,J_R);
j=max([j_L;j_E;j_R]);
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

% 组装各区域的矩阵
[K_L,M_L,BaseIndex_L,DofIndex_L] = AssemblePiezMatWave(Mat1DRef,BaseFunIndex_ref,...
    jLcoef,J_L,j0,MaterialPara,interval_L); % 左
[K_E,M_E,BaseIndex_E,DofIndex_E] = AssemblePiezMatWave(Mat1DRef,BaseFunIndex_ref,...
    jEcoef,J_E,j0,MaterialPara,interval_E); % 中
[K_R,M_R,BaseIndex_R,DofIndex_R] = AssemblePiezMatWave(Mat1DRef,BaseFunIndex_ref,...
    jRcoef,J_R,j0,MaterialPara,interval_R); % 右
% 处理Dirichlet边界条件

%  中
K_E=K_E-omega_p^2*M_E;
BottomIndexE=find(BaseIndex_E(:,18)==interval_E(3,1));
TopIndexE=find(BaseIndex_E(:,18)==interval_E(3,2));
DiriInfo_E={1,BottomIndexE,0,0;2,BottomIndexE,0,0;3,BottomIndexE,0,0;4,BottomIndexE,0,0;...
    4,TopIndexE,@(x,y,z)ones(size(x,1),size(x,2))*V0,1};
[K_E,F_E,uIndexI_E,vIndexI_E,wIndexI_E,phiIndexI_E,uB_E,vB_E,wB_E,phiB_E] ...
    = TreatPiezDiriBoundWave3D(K_E,zeros(size(K_E,1),1),DiriInfo_E,BaseIndex_E,WaveBaseType);
BaseIndexTopE=BaseIndex_E(phiB_E(:,1),16:18); % 中间区域上表面基函数的坐标索引,行的含义与phiB_E一致
%  左
K_L=K_L-omega_p^2*M_L;
%   提取右上脚边
N_L=size(BaseIndex_L,1);
[TRIndexL,TRIndexL2E]=find(BaseIndex_L(:,16)==BaseIndexTopE(:,1)'&...
    BaseIndex_L(:,17)==BaseIndexTopE(:,2)'&...
    BaseIndex_L(:,18)==BaseIndexTopE(:,3)'); % TRIndexL：左区域中的右上角边基函数编号；TRIndexL2E：对应的BaseIndexTopE的行
phiB_L=[TRIndexL,phiB_E(TRIndexL2E,2)]; % 右上脚边的数值解
phiB_L(BaseIndex_L(phiB_L(:,1),18)==interval_L(3,1),:)=[];
BottomIndexL=find(BaseIndex_L(:,18)==interval_L(3,1)); % 底面的基函数编号
F_L=zeros(size(K_L,1),1);
F_L=F_L-sum(phiB_L(:,2).'.*K_L(:,phiB_L(:,1)+3*N_L),2);
%    边界处基函数的数值解，1：对应BaseIndex_L的行；2：数值解
uB_L=[BottomIndexL,zeros(length(BottomIndexL),1)];
vB_L=[BottomIndexL,zeros(length(BottomIndexL),1)];
wB_L=[BottomIndexL,zeros(length(BottomIndexL),1)];
phiB_L=[BottomIndexL,zeros(length(BottomIndexL),1);phiB_L];
%    内部基函数的索引，1：对应刚度矩阵K_L的位置；2：对应BaseIndex_L的行
uIndexI_L=(1:N_L)'; vIndexI_L=uIndexI_L; wIndexI_L=uIndexI_L; phiIndexI_L=uIndexI_L;
uIndexI_L(uB_L(:,1))=[]; 
vIndexI_L(vB_L(:,1))=[]; 
wIndexI_L(wB_L(:,1))=[]; 
phiIndexI_L(phiB_L(:,1))=[];
uIndexI_L=[(1:length(uIndexI_L))',uIndexI_L]; 
vIndexI_L=[size(uIndexI_L,1)+(1:length(vIndexI_L))',vIndexI_L];
wIndexI_L=[size(uIndexI_L,1)+size(vIndexI_L,1)+(1:length(wIndexI_L))',wIndexI_L];
phiIndexI_L=[size(uIndexI_L,1)+size(vIndexI_L,1)+size(wIndexI_L,1)+(1:length(phiIndexI_L))',phiIndexI_L];
index_del_L=[uB_L(:,1);vB_L(:,1)+N_L;wB_L(:,1)+2*N_L;phiB_L(:,1)+3*N_L];
K_L(index_del_L,:)=[]; K_L(:,index_del_L)=[]; F_L(index_del_L)=[];

%  右
K_R=K_R-omega_p^2*M_R;
%   提取左上脚边
N_R=size(BaseIndex_R,1);
[TRIndexR,TRIndexR2E]=find(BaseIndex_R(:,16)==BaseIndexTopE(:,1)'&...
    BaseIndex_R(:,17)==BaseIndexTopE(:,2)'&...
    BaseIndex_R(:,18)==BaseIndexTopE(:,3)'); % TRIndexR：左区域中的左上角边基函数编号；TRIndexR2E：对应的BaseIndexTopE的行
phiB_R=[TRIndexR,phiB_E(TRIndexR2E,2)]; % 右上脚边的数值解
phiB_R(BaseIndex_R(phiB_R(:,1),18)==interval_R(3,1),:)=[];
BottomIndexR=find(BaseIndex_R(:,18)==interval_R(3,1)); % 底面的基函数编号
F_R=zeros(size(K_R,1),1);
F_R=F_R-sum(phiB_R(:,2).'.*K_R(:,phiB_R(:,1)+3*N_R),2);
%    边界处基函数的数值解，1：对应BaseIndex_L的行；2：数值解
uB_R=[BottomIndexR,zeros(length(BottomIndexR),1)];
vB_R=[BottomIndexR,zeros(length(BottomIndexR),1)];
wB_R=[BottomIndexR,zeros(length(BottomIndexR),1)];
phiB_R=[BottomIndexR,zeros(length(BottomIndexR),1);phiB_R];
%    内部基函数的索引，1：对应刚度矩阵K_L的位置；2：对应BaseIndex_L的行
uIndexI_R=(1:N_R)'; vIndexI_R=uIndexI_R; wIndexI_R=uIndexI_R; phiIndexI_R=uIndexI_R;
uIndexI_R(uB_R(:,1))=[]; vIndexI_R(vB_R(:,1))=[]; wIndexI_R(wB_R(:,1))=[]; phiIndexI_R(phiB_R(:,1))=[];
uIndexI_R=[(1:length(uIndexI_R))',uIndexI_R]; 
vIndexI_R=[size(uIndexI_R,1)+(1:length(vIndexI_R))',vIndexI_R];
wIndexI_R=[size(uIndexI_R,1)+size(vIndexI_R,1)+(1:length(wIndexI_R))',wIndexI_R];
phiIndexI_R=[size(uIndexI_R,1)+size(vIndexI_R,1)+size(wIndexI_R,1)+(1:length(phiIndexI_R))',phiIndexI_R];
index_del_R=[uB_R(:,1);vB_R(:,1)+N_R;wB_R(:,1)+2*N_R;phiB_R(:,1)+3*N_R];
K_R(index_del_R,:)=[]; K_R(:,index_del_R)=[]; F_R(index_del_R)=[];

%  左
% K_L=K_L-omega_p^2*M_L;
% BottomIndexL=find(BaseIndex_L(:,18)==interval_L(3,1));
% DiriInfo_L={1,BottomIndexL,0,0;2,BottomIndexL,0,0;3,BottomIndexL,0,0;4,BottomIndexL,0,0};
% [K_L,F_L,uIndexI_L,vIndexI_L,wIndexI_L,phiIndexI_L,uB_L,vB_L,wB_L,phiB_L] ...
%     = TreatPiezDiriBoundWave3D(K_L,zeros(size(K_L,1),1),DiriInfo_L,BaseIndex_L,WaveBaseType);

%  右
% K_R=K_R-omega_p^2*M_R;
% BottomIndexR=find(BaseIndex_R(:,18)==interval_R(3,1));
% DiriInfo_R={1,BottomIndexR,0,0;2,BottomIndexR,0,0;3,BottomIndexR,0,0;4,BottomIndexR,0,0};
% [K_R,F_R,uIndexI_R,vIndexI_R,wIndexI_R,phiIndexI_R,uB_R,vB_R,wB_R,phiB_R] ...
%     = TreatPiezDiriBoundWave3D(K_R,zeros(size(K_R,1),1),DiriInfo_R,BaseIndex_R,WaveBaseType);


% 构造用于生成Lagrange乘子矩阵的信息表
LagInfo_L=genLagInfo(BaseIndex_L,uIndexI_L,vIndexI_L,wIndexI_L,phiIndexI_L,interval_L);
LagInfo_E=genLagInfo(BaseIndex_E,uIndexI_E,vIndexI_E,wIndexI_E,phiIndexI_E,interval_E);
LagInfo_R=genLagInfo(BaseIndex_R,uIndexI_R,vIndexI_R,wIndexI_R,phiIndexI_R,interval_R);
% 组装连接矩阵
% 左、中
[Mat_L,Mat_E_L] = Connect2WaveByDualMortar(LagInfo_L,LagInfo_E);
% 中、右
[Mat_R,Mat_E_R] = Connect2WaveByDualMortar(LagInfo_R,LagInfo_E);
% 拼装矩阵
B=[Mat_L,sparse(size(Mat_L,1),size(Mat_E_R,2));...
    -Mat_E_L,Mat_E_R;...
    sparse(size(Mat_R,1),size(Mat_E_L,2)),-Mat_R];
K=blkdiag(K_L,K_E,K_R);
K=[K,B;B',sparse(size(B,2),size(B,2))];
F=[F_L;F_E;F_R;zeros(size(B,2),1)];
% 计算数值解
uh1=K\F;
% 生成各子区域的数值解
%  左
[u_L,v_L,w_L,phi_L]=genNumSol(uh1,uIndexI_L,uB_L,vIndexI_L,vB_L,wIndexI_L,wB_L,phiIndexI_L,phiB_L);
u_L=@(x,y,z) ApproxWaveFun3D(x,y,z,u_L,BaseIndex_L,[0,0,0],WaveBaseType);
v_L=@(x,y,z) ApproxWaveFun3D(x,y,z,v_L,BaseIndex_L,[0,0,0],WaveBaseType);
w_L=@(x,y,z) ApproxWaveFun3D(x,y,z,w_L,BaseIndex_L,[0,0,0],WaveBaseType);
phi_L=@(x,y,z) ApproxWaveFun3D(x,y,z,phi_L,BaseIndex_L,[0,0,0],WaveBaseType);
%  中
[u_E,v_E,w_E,phi_E]=genNumSol(uh1,uIndexI_E,uB_E,vIndexI_E,vB_E,wIndexI_E,wB_E,phiIndexI_E,phiB_E);
u_E=@(x,y,z) ApproxWaveFun3D(x,y,z,u_E,BaseIndex_E,[0,0,0],WaveBaseType);
v_E=@(x,y,z) ApproxWaveFun3D(x,y,z,v_E,BaseIndex_E,[0,0,0],WaveBaseType);
w_E=@(x,y,z) ApproxWaveFun3D(x,y,z,w_E,BaseIndex_E,[0,0,0],WaveBaseType);
phi_E=@(x,y,z) ApproxWaveFun3D(x,y,z,phi_E,BaseIndex_E,[0,0,0],WaveBaseType);
%  右
[u_R,v_R,w_R,phi_R]=genNumSol(uh1,uIndexI_R,uB_R,vIndexI_R,vB_R,wIndexI_R,wB_R,phiIndexI_R,phiB_R);
u_R=@(x,y,z) ApproxWaveFun3D(x,y,z,u_R,BaseIndex_R,[0,0,0],WaveBaseType);
v_R=@(x,y,z) ApproxWaveFun3D(x,y,z,v_R,BaseIndex_R,[0,0,0],WaveBaseType);
w_R=@(x,y,z) ApproxWaveFun3D(x,y,z,w_R,BaseIndex_R,[0,0,0],WaveBaseType);
phi_R=@(x,y,z) ApproxWaveFun3D(x,y,z,phi_R,BaseIndex_R,[0,0,0],WaveBaseType);
% 绘制位移图
figure(1)
DrawCubeDomain(interval_L);
DrawCubeDomain(interval_E);
DrawCubeDomain(interval_R);
hold on
% 左
stepx=(interval_L(1,2)-interval_L(1,1))/20;
stepy=(interval_L(2,2)-interval_L(2,1))/10;
stepz=(interval_L(3,2)-interval_L(3,1))/100;
[xx,yy,zz]=genDiscretePoints(interval_L,stepx,stepy,stepz);
u_L=u_L(xx,yy,zz); v_L=v_L(xx,yy,zz); w_L=w_L(xx,yy,zz);
Pos_L=[xx,yy,zz]+[u_L,v_L,w_L]*5e3;
Disp_L=sqrt(u_L.^2+v_L.^2+w_L.^2);
scatter3(Pos_L(:,1),Pos_L(:,2),Pos_L(:,3),10,Disp_L,'filled')
% 中
stepx=(interval_E(1,2)-interval_E(1,1))/20;
stepy=(interval_E(2,2)-interval_E(2,1))/10;
stepz=(interval_E(3,2)-interval_E(3,1))/100;
[xx,yy,zz]=genDiscretePoints(interval_E,stepx,stepy,stepz);
u_E=u_E(xx,yy,zz); v_E=v_E(xx,yy,zz); w_E=w_E(xx,yy,zz);
Pos_E=[xx,yy,zz]+[u_E,v_E,w_E]*5e3;
Disp_E=sqrt(u_E.^2+v_E.^2+w_E.^2);
scatter3(Pos_E(:,1),Pos_E(:,2),Pos_E(:,3),10,Disp_E,'filled')
% 右
stepx=(interval_R(1,2)-interval_R(1,1))/20;
stepy=(interval_R(2,2)-interval_R(2,1))/10;
stepz=(interval_R(3,2)-interval_R(3,1))/100;
[xx,yy,zz]=genDiscretePoints(interval_R,stepx,stepy,stepz);
u_R=u_R(xx,yy,zz); v_R=v_R(xx,yy,zz); w_R=w_R(xx,yy,zz);
Pos_R=[xx,yy,zz]+[u_R,v_R,w_R]*5e3;
Disp_R=sqrt(u_R.^2+v_R.^2+w_R.^2);
scatter3(Pos_R(:,1),Pos_R(:,2),Pos_R(:,3),10,Disp_R,'filled')
colormap jet

% 绘制电势图
figure(2)
DrawCubeDomain(interval_L);
DrawCubeDomain(interval_E);
DrawCubeDomain(interval_R);
hold on
% 左
stepx=(interval_L(1,2)-interval_L(1,1))/20;
stepy=(interval_L(2,2)-interval_L(2,1))/10;
stepz=(interval_L(3,2)-interval_L(3,1))/100;
[xx,yy,zz]=genDiscretePoints(interval_L,stepx,stepy,stepz);
phi_L=phi_L(xx,yy,zz);
scatter3(xx,yy,zz,10,phi_L,'filled')
% 中
stepx=(interval_E(1,2)-interval_E(1,1))/20;
stepy=(interval_E(2,2)-interval_E(2,1))/10;
stepz=(interval_E(3,2)-interval_E(3,1))/100;
[xx,yy,zz]=genDiscretePoints(interval_E,stepx,stepy,stepz);
phi_E=phi_E(xx,yy,zz);
scatter3(xx,yy,zz,10,phi_E,'filled')
% 右
stepx=(interval_R(1,2)-interval_R(1,1))/20;
stepy=(interval_R(2,2)-interval_R(2,1))/10;
stepz=(interval_R(3,2)-interval_R(3,1))/100;
[xx,yy,zz]=genDiscretePoints(interval_R,stepx,stepy,stepz);
phi_R=phi_R(xx,yy,zz);
scatter3(xx,yy,zz,10,phi_R,'filled')
colormap jet

function LagInfo=genLagInfo(BaseIndex,uIndexI,vIndexI,wIndexI,phiIndexI,interval)
BaseIndex_u=BaseIndex(uIndexI(:,2),:);
BaseIndex_v=BaseIndex(vIndexI(:,2),:);
BaseIndex_w=BaseIndex(wIndexI(:,2),:);
BaseIndex_phi=BaseIndex(phiIndexI(:,2),:);
N_u=size(BaseIndex_u,1);
N_v=size(BaseIndex_v,1);
N_w=size(BaseIndex_w,1);
N_phi=size(BaseIndex_phi,1);
LagInfo={interval,...
    [BaseIndex_u(:,4:6);BaseIndex_v(:,4:6);BaseIndex_w(:,4:6);BaseIndex_phi(:,4:6)],...
    [BaseIndex_u(:,7:9);BaseIndex_v(:,7:9);BaseIndex_w(:,7:9);BaseIndex_phi(:,7:9)],...
    [BaseIndex_u(:,1:3);BaseIndex_v(:,1:3);BaseIndex_w(:,1:3);BaseIndex_phi(:,1:3)],...
    [BaseIndex_u(:,16:18);BaseIndex_v(:,16:18);BaseIndex_w(:,16:18);BaseIndex_phi(:,16:18)],...
    [ones(N_u,1);2*ones(N_v,1);3*ones(N_w,1);4*ones(N_phi,1)]};
end
function [u,v,w,phi]=genNumSol(uh,uIndexI,uB,vIndexI,vB,wIndexI,wB,phiIndexI,phiB)
u=sortrows([uIndexI(:,2),uh(uIndexI(:,1));uB]);
v=sortrows([vIndexI(:,2),uh(vIndexI(:,1));vB]);
w=sortrows([wIndexI(:,2),uh(wIndexI(:,1));wB]);
phi=sortrows([phiIndexI(:,2),uh(phiIndexI(:,1));phiB]);
u=u(:,2);
v=v(:,2);
w=w(:,2);
phi=phi(:,2);
end
function [xx,yy,zz]=genDiscretePoints(interval,stepx,stepy,stepz)
xx=interval(1,1):stepx:interval(1,2);
yy=interval(2,1):stepy:interval(2,2);
zz=interval(3,1):stepz:interval(3,2);
[xx,yy,zz]=meshgrid(xx,yy,zz);
xx=xx(:);
yy=yy(:);
zz=zz(:);
end
