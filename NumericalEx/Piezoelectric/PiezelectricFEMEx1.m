% 具有Dirichlet边界条件的算例
clear
% clc
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
addpath(genpath(currentPath));
% 材料参数
ModelCoeff = 'D:\Code\M\Mortar_FEM_Wavelet\Piezoelectric\Data\ModelCoefExample.mat';
%% 构造方程
syms u(x,y,z) v(x,y,z) w(x,y,z) phi(x,y,z)
% EX9(J=9时有超收敛)
u(x,y,z) = sin(x + y + z)+1;
v(x,y,z) = sin(pi*x).*cos(pi*y).*sin(pi*z)+1;
w(x,y,z) = x.*sin(y)+y.*sin(z)+z.*sin(x);
phi(x,y,z) = sin(x+y.^2+z.^3)+1;
u_exact = matlabFunction(u);
v_exact = matlabFunction(v);
w_exact = matlabFunction(w);
phi_exact = matlabFunction(phi);
[f1,f2,f3,g] = GenPiezoeleModel(ModelCoeff,u,v,w,x,y,z,phi);
f1 =@(x,y,z)f1(x,y,z).*ones(size(x,1),1);
f2 =@(x,y,z)f2(x,y,z).*ones(size(x,1),1);
f3 =@(x,y,z)f3(x,y,z).*ones(size(x,1),1);
g =@(x,y,z) g(x,y,z).*ones(size(x,1),1);
%% 数值求解
% 导入材料参数
load(ModelCoeff)
c_LN = cell2mat(materials(2));
e_LN = cell2mat(materials(5));
epcl0 = cell2mat(materials(6));
epcl_LN = cell2mat(materials(7));
epcl_LN = epcl_LN*epcl0; 
% 几何参数
ori1 = [0,0,0];
a1 = 1;b1=1;c1=1;
% 网格参数
% type="quadratic";
type="linear";
N = 4;
Nx1=N+1;Ny1=N+1;Nz1=N+1;
% 生成压电方程矩阵与载荷矢量
[K,~,F,Dof_Index] =...
    AssemblePiezMatFEM(c_LN,e_LN,epcl_LN,0,...
    ori1,a1,b1,c1,Nx1,Ny1,Nz1,...
    [],[],[],type,...
    {f1,f2,f3,g});
% 处理Dirichlete边界条件
boundary_u=find(Dof_Index(:,1)==1&...
    (Dof_Index(:,2)==ori1(1)|Dof_Index(:,2)==ori1(1)+a1|...
    Dof_Index(:,3)==ori1(2)|Dof_Index(:,3)==ori1(2)+b1|...
    Dof_Index(:,4)==ori1(3)|Dof_Index(:,4)==ori1(3)+c1));
boundary_v=find(Dof_Index(:,1)==2&...
    (Dof_Index(:,2)==ori1(1)|Dof_Index(:,2)==ori1(1)+a1|...
    Dof_Index(:,3)==ori1(2)|Dof_Index(:,3)==ori1(2)+b1|...
    Dof_Index(:,4)==ori1(3)|Dof_Index(:,4)==ori1(3)+c1));
boundary_w=find(Dof_Index(:,1)==3&...
    (Dof_Index(:,2)==ori1(1)|Dof_Index(:,2)==ori1(1)+a1|...
    Dof_Index(:,3)==ori1(2)|Dof_Index(:,3)==ori1(2)+b1|...
    Dof_Index(:,4)==ori1(3)|Dof_Index(:,4)==ori1(3)+c1));
boundary_phi=find(Dof_Index(:,1)==4&...
    (Dof_Index(:,2)==ori1(1)|Dof_Index(:,2)==ori1(1)+a1|...
    Dof_Index(:,3)==ori1(2)|Dof_Index(:,3)==ori1(2)+b1|...
    Dof_Index(:,4)==ori1(3)|Dof_Index(:,4)==ori1(3)+c1));

F=F-sum(u_exact(Dof_Index(boundary_u,2),Dof_Index(boundary_u,3),Dof_Index(boundary_u,4)).'.*K(:,boundary_u),2)-...
    sum(v_exact(Dof_Index(boundary_v,2),Dof_Index(boundary_v,3),Dof_Index(boundary_v,4)).'.*K(:,boundary_v),2)-...
    sum(w_exact(Dof_Index(boundary_w,2),Dof_Index(boundary_w,3),Dof_Index(boundary_w,4)).'.*K(:,boundary_w),2)-...
    sum(phi_exact(Dof_Index(boundary_phi,2),Dof_Index(boundary_phi,3),Dof_Index(boundary_phi,4)).'.*K(:,boundary_phi),2);
K([boundary_u;boundary_v;boundary_w;boundary_phi],:) = [];
K(:,[boundary_u;boundary_v;boundary_w;boundary_phi]) = [];
F([boundary_u;boundary_v;boundary_w;boundary_phi],:) = [];
Dof_Index([boundary_u;boundary_v;boundary_w;boundary_phi],:)=[];
% 解方程
SolNum = K\F;
SolNum = reshape(SolNum,[],4);
disp('相对误差')
Pb_sub=Dof_Index(Dof_Index(:,1)==1,2:4);
Solexact = ...
    [u_exact(Pb_sub(:,1),Pb_sub(:,2),Pb_sub(:,3)),...
    v_exact(Pb_sub(:,1),Pb_sub(:,2),Pb_sub(:,3)),...
    w_exact(Pb_sub(:,1),Pb_sub(:,2),Pb_sub(:,3)),...
    phi_exact(Pb_sub(:,1),Pb_sub(:,2),Pb_sub(:,3))];
err_all = norm(SolNum-Solexact,inf)./norm(Solexact,inf)
