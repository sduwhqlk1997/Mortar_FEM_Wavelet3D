% 用HCT方法求解交界面电荷
clear
% clc
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
% 材料参数
ModelCoeff = 'D:\Code\M\Mortar_FEM_Wavelet\Piezoelectric\Data\ModelCoef2.mat';
load(ModelCoeff,'materials','Geo')
omega = 4.39822971502571*1e9;
c_LN = cell2mat(materials(2));
lamda_Al = cell2mat(materials(3));
mu_Al = cell2mat(materials(4));
e_LN = cell2mat(materials(5));
epcl_LN = cell2mat(materials(7))*cell2mat(materials(6));
density_LN = cell2mat(materials(8));
density_Al = cell2mat(materials(9));
% 级联参数
N_HCT=5; % 级联次数
N_IDT=2^N_HCT;
% 几何参数
ori1 = cell2mat(Geo(1));a1 = cell2mat(Geo(2));b1 = cell2mat(Geo(3));c1 = cell2mat(Geo(4));
ori2 = cell2mat(Geo(5));a2 = cell2mat(Geo(6));b2 = cell2mat(Geo(7));c2 = cell2mat(Geo(8));
interval_sub = [ori1.',ori1.'+[a1;b1;c1]];
interval_ele = [ori2.',ori2.'+[a2;b2;c2]];
% PML设置
% d_max=1;
d_max=1e6/omega;
d_pml = 2e-6;
n_pml=2;
x_p_L = ori1(1)-d_pml; x_p_R = ori1(1)+a1*N_IDT+d_pml;
x_p_B = ori1(3)-d_pml;
Df_L = @(x) d_max*(1-(x-x_p_L).^2/d_pml^2).^n_pml; % 左PML
Df_R = @(x) d_max*(1-(x-x_p_R).^2/d_pml^2).^n_pml; % 右PML
Df_B = @(x) d_max*(1-(x-x_p_B).^2/d_pml^2).^n_pml; % 下PML
% Df_L = @(x) zeros(size(x,1),size(x,2)); % 左PML
% Df_R = @(x) zeros(size(x,1),size(x,2)); % 右PML
% Df_B = @(x) zeros(size(x,1),size(x,2)); % 下PML
% 网格参数
N=8;
Nx1 = 2*N+1; Ny1 = 2; Nz1 = 2*N+1;
Nx2 = N+1; Nz2 = N+1;
Nx_pml=N+1;

% Nx1 = 17; Ny1 = 2; Nz1 = 17;
% Nx2 = 9; Nz2 = 5;
% Nx_pml=5;

%%%%%%%%%%%测试用%%%%%%%%%%%%
% Nx1 = 9; Ny1 = 2; Nz1 = 9;
% Nx2 = 5; Nz2 = 5;
% Nx_pml=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type="quadratic";
% type="linear";
% 组装单根IDT刚度矩阵
[K_IDT,M_IDT,Dof_Index_IDT,~,~,~,xStep_sub,yStep_sub,zStep_sub] = AssembleSAWMatFEM(c_LN,e_LN,epcl_LN,density_LN,...
    lamda_Al,mu_Al,density_Al,...
    ori1,a1,b1,c1,Nx1,Ny1,Nz1,...
    ori2,a2,b2,c2,Nx2,Nz2,true,type);
K_IDT=K_IDT-omega^2*M_IDT;
% 组装PML区域矩阵
[K_L,M_L,Dof_index_L] = AssemblePiezPMLMatFEM(Df_L,@(y) zeros(size(y,1),1),@(z) zeros(size(z,1),1),...
    c_LN,e_LN,epcl_LN,density_LN,...
    [x_p_L,ori1(2:3)],d_pml,b1,c1,Nx_pml,Ny1,Nz1,...
    [],yStep_sub,zStep_sub,type);
K_L=K_L-omega^2*M_L;

[K_R,M_R,Dof_index_R] = AssemblePiezPMLMatFEM(Df_R,@(y) zeros(size(y,1),1),@(z) zeros(size(z,1),1),...
    c_LN,e_LN,epcl_LN,density_LN,...
    [x_p_R-d_pml,ori1(2:3)],d_pml,b1,c1,Nx_pml,Ny1,Nz1,...
    [],yStep_sub,zStep_sub,type);
K_R=K_R-omega^2*M_R;

[K_B,M_B,Dof_index_B] = AssemblePiezPMLMatFEM(@(x) zeros(size(x,1),1),@(y) zeros(size(y,1),1),Df_B,...
    c_LN,e_LN,epcl_LN,density_LN,...
    [ori1(1:2),ori1(3)-d_pml],a1,b1,d_pml,Nx1,Ny1,Nx_pml,...
    xStep_sub,yStep_sub,[],type);
K_B=K_B-omega^2*M_B;

[K_LB,M_LB,Dof_index_LB] = AssemblePiezPMLMatFEM(Df_L,@(y) zeros(size(y,1),1),Df_B,...
    c_LN,e_LN,epcl_LN,density_LN,...
    [x_p_L,ori1(2),ori1(3)-d_pml],d_pml,b1,d_pml,Nx_pml,Ny1,Nx_pml,...
    [],yStep_sub,[],type);
K_LB=K_LB-omega^2*M_LB;

[K_RB,M_RB,Dof_index_RB] = AssemblePiezPMLMatFEM(Df_R,@(y) zeros(size(y,1),1),Df_B,...
    c_LN,e_LN,epcl_LN,density_LN,...
    [x_p_R-d_pml,ori1(2),ori1(3)-d_pml],d_pml,b1,d_pml,Nx_pml,Ny1,Nx_pml,...
    [],yStep_sub,[],type);
K_RB=K_RB-omega^2*M_RB;
% 拼接矩阵
[K_IDT,~,Dof_Index_IDT] = Connect2RectRegion(K_IDT,K_B,Dof_Index_IDT,Dof_index_B);
interval_sub(3,1)=x_p_B;
[K_L,~,Dof_index_L] = Connect2RectRegion(K_L,K_LB,Dof_index_L,Dof_index_LB);
interval_L=[x_p_L,ori1(1);interval_sub(2,:);x_p_B,interval_sub(3,2)];
[K_R,~,Dof_index_R] = Connect2RectRegion(K_R,K_RB,Dof_index_R,Dof_index_RB);
interval_R=[x_p_R-d_pml,x_p_R;interval_sub(2,:);x_p_B,interval_sub(3,2)];
% 处理周期边界条件
[K_IDT,Dof_Index_IDT]=TreatPeriodBound(K_IDT,Dof_Index_IDT,interval_sub);
[K_L,Dof_index_L]=TreatPeriodBound(K_L,Dof_index_L,interval_L);
[K_R,Dof_index_R]=TreatPeriodBound(K_R,Dof_index_R,interval_R);
% 处理边界条件
%  IDT
bottom=abs(Dof_Index_IDT(:,4)-x_p_B)/c1<1e-10;
K_IDT(bottom,:)=[];
K_IDT(:,bottom)=[];
Dof_Index_IDT(bottom,:)=[];
%  左PML
bottom=abs(Dof_index_L(:,2)-x_p_L)/d_pml<1e-10|abs(Dof_index_L(:,4)-x_p_B)/d_pml<1e-10;
K_L(bottom,:)=[];
K_L(:,bottom)=[];
Dof_index_L(bottom,:)=[];
%  右PML
bottom=abs(Dof_index_R(:,2)-x_p_R)/d_pml<1e-10|abs(Dof_index_R(:,4)-x_p_B)/d_pml<1e-10;
K_R(bottom,:)=[];
K_R(:,bottom)=[];
Dof_index_R(bottom,:)=[];
% 计算schur补

[K_IDT,Dof_Index_IDT] = IDT_Schur(K_IDT,Dof_Index_IDT,interval_sub,interval_ele,true);
[K_L,K_R,Dof_index_L,Dof_index_R] =...
    PML_Schur(K_L,K_R,Dof_index_L,Dof_index_R,interval_L,interval_R);
% HCT
[K_IDT,Dof_Index_IDT] = IDT_HCT(K_IDT,Dof_Index_IDT,N_HCT);
[K_IDT,Dof_Index_IDT] ...
    = CombineIDTandPML(K_IDT,Dof_Index_IDT,K_L,K_R,Dof_index_L,Dof_index_R);
% 计算表面电荷
V=1;
index_L=Dof_Index_IDT(:,1)==-1;
index_R=Dof_Index_IDT(:,1)==-2;
index_b=index_L|index_R;
index_e = Dof_Index_IDT(:,1)>0;
Q=(K_IDT(index_e,index_b)*(K_IDT(index_b,index_b)\sum(K_IDT(index_b,index_e),2)*V)-...
    sum(K_IDT(index_e,index_e),2)*V)*400;
Ad=omega*sum(Q)*1j
% Q=K_IDT(index_e,index_e)-K_IDT(index_e,index_b)*(K_IDT(index_b,index_b)\K_IDT(index_b,index_e))
function [K,DofIndex]=TreatPeriodBound(K,DofIndex,interval)
index_front=find(DofIndex(:,3)==interval(2,1));
index_back=find(DofIndex(:,3)==interval(2,2));
[i,j]=find(DofIndex(index_front,1)==DofIndex(index_back,1).'&...
    abs(DofIndex(index_front,2)-DofIndex(index_back,2).')/(interval(1,2)-interval(1,1))<1e-10&...
    abs(DofIndex(index_front,4)-DofIndex(index_back,4).')/(interval(3,2)-interval(3,1))<1e-10);% i 前面，j 后面
i=index_front(i);
j=index_back(j);
K(i,:)=K(i,:)+K(j,:);
K(:,i)=K(:,i)+K(:,j);
K(j,:)=[];
K(:,j)=[];
DofIndex(j,:)=[];
end