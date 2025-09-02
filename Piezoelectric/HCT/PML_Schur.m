function [K_L,K_R,Dof_indexHCT_L,Dof_indexHCT_R,SizeSchurPML] =...
    PML_Schur(K_L,K_R,Dof_index_L,Dof_index_R,interval_L,interval_R)
%PML_SCHUR 计算左右PML的Schur补，并生成用于HCT的索引矩阵
%   K_L,K_R：左、右PML的刚度矩阵
%   Dof_index_L,Dof_index_R：左、右PML的索引
%       Dof_index：有限元索引矩阵
%           第一列自由度类型1:u,2:v,3:w,4:phi；2-4列坐标
%   Dof_indexHCT：HCT索引矩阵
%       第一列：-1：左边界，-2：右边界，1,2,3,....：电极电自由度，对应电极的编号
%       第二、三列：y、z坐标
%       第四列：自由度类型
%   注：在使用本函数之前PML区域的矩阵应处理完所有边界条件

index_b_L=abs(Dof_index_L(:,2)-interval_L(1,2))/(interval_L(1,2)-interval_L(1,1))<1e-10; % 左PML的右边界索引
index_i_L=~index_b_L; % 左PML的内部自由度索引

index_b_R=abs(Dof_index_R(:,2)-interval_R(1,1))/(interval_R(1,2)-interval_R(1,1))<1e-10; % 右PML的左边界索引
index_i_R=~index_b_R; % 右PML的内部自由度索引

% 生成HCT索引
Dof_indexHCT_L = [zeros(size(Dof_index_L,1),1),Dof_index_L(:,3),Dof_index_L(:,4),Dof_index_L(:,1)];
Dof_indexHCT_L(index_b_L,1)=-2;

Dof_indexHCT_R = [zeros(size(Dof_index_R,1),1),Dof_index_R(:,3),Dof_index_R(:,4),Dof_index_R(:,1)];
Dof_indexHCT_R(index_b_R,1)=-1;

% 消去内部自由度
SizeSchurPML=size(K_L(index_i_L,index_i_L),1);
Kbb=K_L(index_i_L,index_i_L)\K_L(index_i_L,index_b_L);
Kbb=K_L(index_b_L,index_i_L)*Kbb;
K_L=K_L(index_b_L,index_b_L)-Kbb;
% Dof_indexHCT_L(index_i_L,:)=[];
Dof_indexHCT_L=Dof_indexHCT_L(index_b_L,:);

SizeSchurPML=[SizeSchurPML,size(K_R(index_i_R,index_i_R),1)];
Kbb=K_R(index_i_R,index_i_R)\K_R(index_i_R,index_b_R);
Kbb=K_R(index_b_R,index_i_R)*Kbb;
K_R=K_R(index_b_R,index_b_R)-Kbb;
Dof_indexHCT_R=Dof_indexHCT_R(index_b_R,:);
end

