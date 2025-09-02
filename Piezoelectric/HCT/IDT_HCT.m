function [K,Dof_indexHCT] = IDT_HCT(K,Dof_indexHCT,N_iter)
%IDT_HCT 级联所有IDT（不含左、右PML）
% 输出
%   Dof_indexHCT级联操作索引矩阵，每行与K每行对应
%       第一列：-1：左边界，-2：右边界，1,2,3,....：电极电自由度，对应电极的编号
%       第二、三列：y、z坐标
%       第三列：自由度类型
% 输入
%   K：单根IDT的Schur补矩阵
%   N_iter：级联次数

index_L=Dof_indexHCT(:,1)==-1;
index_R=Dof_indexHCT(:,1)==-2;
index_e=~(index_L|index_R);
for i=1:N_iter
    N=size(K,1);
    K=blkdiag(K,K);
    index_L1=[index_L;false(N,1)];
    index_R1=[index_R;false(N,1)];
    index_e1=[index_e;false(N,1)];
    index_L2=[false(N,1);index_L];
    index_R2=[false(N,1);index_R];
    index_e2=[false(N,1);index_e];
    Dof_indexHCT = repmat(Dof_indexHCT,2,1);
    K(index_R1,:) = K(index_R1,:)+K(index_L2,:);
    K(:,index_R1) = K(:,index_R1)+K(:,index_L2);
    K(index_L2,:)=[];
    K(:,index_L2)=[];
    index_L1(index_L2,:)=[];
    index_R1(index_L2,:)=[];
    index_R2(index_L2,:)=[];
    index_e1(index_L2,:)=[];
    index_e2(index_L2,:)=[];
    Dof_indexHCT(index_L2,:)=[];
    index_L=index_L1;
    index_R=index_R2;
    index_b=index_L|index_R|index_e1|index_e2;
    index_i=index_R1;
    Kbb=K(index_i,index_i)\K(index_i,index_b);
    Kbb=K(index_b,index_i)*Kbb;
    K=K(index_b,index_b)-Kbb;
    % 更新自由度索引
    Dof_indexHCT=Dof_indexHCT(index_b,:);
    index_L=index_L(index_b,:);
    index_R=index_R(index_b,:);
    index_e1=index_e1(index_b,:);
    index_e2=index_e2(index_b,:);
    % Dof_indexHCT(index_e2,1)=Dof_indexHCT(index_e2,1)+size(Dof_indexHCT(index_e1,1),1);
    Dof_indexHCT(index_e2,1)=Dof_indexHCT(index_e2,1)+2^(i-1);
    index_e = index_e1|index_e2;
end
end

