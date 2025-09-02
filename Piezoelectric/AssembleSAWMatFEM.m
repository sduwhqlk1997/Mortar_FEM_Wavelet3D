function [K,M,Dof_Index,K_e,M_e,Dof_Index_e,xStep_sub,yStep_sub,zStep_sub] = AssembleSAWMatFEM(c_LN,e_LN,epcl_LN,density_LN,...
    lamda_Al,mu_Al,density_Al,...
    ori1,a1,b1,c1,Nx1,Ny1,Nz1,...
    ori2,a2,b2,c2,Nx2,Nz2,Connect,type)
%ASSEMBLESAWMATFEM 生成SAW方程的刚度矩阵K、质量矩阵M
%   Connect: 是否要拼接电极和压电基底
if ~exist('Connect','var') % 默认电极和基底拼接
    Connect=true;
end
if ori2(3)~=ori1(3)+c1||ori2(1)<ori1(1)||ori2(1)+a2>ori1(1)+a1||ori2(2)<ori1(2)||ori2(2)+b2>ori1(2)+b1%%此时电极位置不正确
    msg = '电极位置不正确';
    error(msg)
end
% 生成网格步长向量
[xStep_sub,yStep_sub,zStep_sub,xStep_ele,yStep_ele,zStep_ele] = ...
    genMeshStepOfSAW(ori1,a1,b1,c1,Nx1,Ny1,Nz1,ori2,a2,c2,Nx2,Nz2);
% 组装压电基底区域矩阵
[K_sub,M_sub,~,Dof_Index_sub] =...
    AssemblePiezMatFEM(c_LN,e_LN,epcl_LN,density_LN,...
    ori1,a1,b1,c1,Nx1,Ny1,Nz1,...
    xStep_sub,yStep_sub,zStep_sub,type);
% 组装电极区域矩阵
[K_e,M_e,Dof_Index_e] = AssembleEleMatFEM(lamda_Al,mu_Al,density_Al,...
    ori2,a2,b2,c2,Nx2,Ny1,Nz2,...
    xStep_ele,yStep_ele,zStep_ele,type);
if Connect==true
    [K,~,Dof_Index] = Connect2RectRegion(K_sub,K_e,Dof_Index_sub,Dof_Index_e);
    [M,~,~] = Connect2RectRegion(M_sub,M_e,Dof_Index_sub,Dof_Index_e);
else
    K=K_sub;
    M=M_sub;
    Dof_Index=Dof_Index_sub;
end
end

