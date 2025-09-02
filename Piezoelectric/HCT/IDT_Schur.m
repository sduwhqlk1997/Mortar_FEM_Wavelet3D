function [K,Dof_indexHCT,SizeSchurIDT] = IDT_Schur(K,Dof_index,interval_sub,interval_ele,combinePhi)
%IDT_SCHUR 计算单根IDT的Schur补，并生成用于HCT的索引矩阵
%   Dof_index：有限元索引矩阵
%       第一列自由度类型1:u,2:v,3:w,4:phi；2-4列坐标
%   Dof_indexHCT：HCT索引矩阵
%       第一列：-1：左边界，-2：右边界，1,2,3,....：电极电自由度，对应电极的编号
%       第二、三列：y、z坐标
%       第三列：自由度类型
%   combinePhi：是否合并交界面处电自由度，默认为合并

if ~exist('combinePhi','var')
    combinePhi=true;
end
index_e = (Dof_index(:,2)-interval_ele(1,1))/(interval_ele(1,2)-interval_ele(1,1))>-1e-10 &...
    (Dof_index(:,2)-interval_ele(1,2))/(interval_ele(1,2)-interval_ele(1,1))<1e-10 &...
    abs(Dof_index(:,4)-interval_ele(3,1))/(interval_ele(3,2)-interval_ele(3,1))<1e-10 &...
    Dof_index(:,1)==4;
if combinePhi==true % 将电极与基底交界处的电自由度缩为1个自由度
    index_ee = find(index_e);
    K(index_ee(1),:)=sum(K(index_ee,:),1);
    K(:,index_ee(1))=sum(K(:,index_ee),2);
    K(index_ee(2:end),:)=[];
    K(:,index_ee(2:end))=[];
    index_e(index_ee(2:end),:)=[];
    Dof_index(index_ee(2:end),:)=[];
end
% 提取边界自由度
index_L = abs(Dof_index(:,2)-interval_sub(1,1))/(interval_sub(1,2)-interval_sub(1,1))<1e-10;
index_R = abs(Dof_index(:,2)-interval_sub(1,2))/(interval_sub(1,2)-interval_sub(1,1))<1e-10;
index_b = index_L|index_R|index_e;
index_i = ~index_b;
% 生成HCT索引
Dof_indexHCT = [zeros(size(Dof_index,1),1),Dof_index(:,3),Dof_index(:,4),Dof_index(:,1)];
Dof_indexHCT(index_L,1)=-1;
Dof_indexHCT(index_R,1)=-2;
Dof_indexHCT(index_e,1)=1;
% 消去内部自由度
SizeSchurIDT=size(K(index_i,index_i),1);
Kbb=K(index_i,index_i)\K(index_i,index_b);
Kbb=K(index_b,index_i)*Kbb;
K=K(index_b,index_b)-Kbb;
% Dof_indexHCT(index_i,:)=[];
Dof_indexHCT=Dof_indexHCT(index_b,:);
end

