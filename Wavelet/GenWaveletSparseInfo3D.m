function [DofIndex_sparse,Table_sparse] = GenWaveletSparseInfo3D(j0,J,...
    BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z,jcoef)
%GENWAVELETSPARSEINFO3D 生成3D稀疏基相关信息
%   BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z:分别为x,y,z方向的一维参考基的基函数信息，
%       1：SpaceType，0表示尺度基，1表示小波基
%       2：level
%       3：基函数编号
%       4：支集左端点
%       5：支集右端点
%       6：基函数坐标索引
%       7：基函数编号
%    如果只输入了BaseFunIndex_x, 则默认其他两个方向与其相同
%   j1coef*j1+j2coef*j2+j3coef*j3<=J
%   DofIndex_sparse:稀疏基索引，1~3列分别为对应x,y,z方向一维小波基的编号，与BaseFunIndex的行对应
%   Table_sparse 稀疏基刚度矩阵的非0元素索引

if (~exist("BaseFunIndex_y","var")) && (~exist("BaseFunIndex_z","var"))
    BaseFunIndex_y = BaseFunIndex_x; 
    BaseFunIndex_z = BaseFunIndex_x; 
end
if ~exist("jcoef","var")
    jcoef=[1;1;1];
end
BaseFunIndex_x(BaseFunIndex_x(:,1)==1 | BaseFunIndex_x(:,2)>j0,2) = ...
    BaseFunIndex_x(BaseFunIndex_x(:,1)==1 | BaseFunIndex_x(:,2)>j0,2)+1;
BaseFunIndex_y(BaseFunIndex_y(:,1)==1 | BaseFunIndex_y(:,2)>j0,2) = ...
    BaseFunIndex_y(BaseFunIndex_y(:,1)==1 | BaseFunIndex_y(:,2)>j0,2)+1;
BaseFunIndex_z(BaseFunIndex_z(:,1)==1 | BaseFunIndex_z(:,2)>j0,2) = ...
    BaseFunIndex_z(BaseFunIndex_z(:,1)==1 | BaseFunIndex_z(:,2)>j0,2)+1;

j_min_x = min(BaseFunIndex_x(:,2)); j_max_x = max(BaseFunIndex_x(:,2)); 
j_min_y = min(BaseFunIndex_y(:,2)); j_max_y = max(BaseFunIndex_y(:,2));
j_min_z = min(BaseFunIndex_z(:,2)); j_max_z = max(BaseFunIndex_z(:,2));

J_num_x = sum(BaseFunIndex_x(:,2)==(j_min_x:j_max_x),1); % x方向一维基每个level基函数个数
J_num_y = sum(BaseFunIndex_y(:,2)==(j_min_y:j_max_y),1); % y方向一维基每个level基函数个数
J_num_z = sum(BaseFunIndex_z(:,2)==(j_min_z:j_max_z),1); % z方向一维基每个level基函数个数

[index1,index2,index3]=meshgrid(j_min_x:j_max_x,j_min_y:j_max_y,j_min_z:j_max_z);
index = [index1(:),index2(:),index3(:)];
index(index(:,1)*jcoef(1)+index(:,2)*jcoef(2)+index(:,3)*jcoef(3)>J,:) = []; % 所有被选中的小波空间组合(稀疏空间)
% index(index(:,1)>J|index(:,2)>J|index(:,3)>J,:) = []; % 所有被选中的小波空间组合(全空间)
Dof_num = sum(J_num_x(index(:,1)-j_min_x+1)...
    .*J_num_y(index(:,2)-j_min_y+1)...
    .*J_num_z(index(:,3)-j_min_z+1)); % 总自由度数
DofIndex_sparse = zeros(Dof_num,3);
pos = 1;

for i = 1:size(index,1)
    [index1,index2,index3] = ...
        meshgrid(BaseFunIndex_x(BaseFunIndex_x(:,2)==index(i,1),end),...
        BaseFunIndex_y(BaseFunIndex_y(:,2)==index(i,2),end),...
        BaseFunIndex_z(BaseFunIndex_z(:,2)==index(i,3),end));
    DofIndex_sparse(pos:pos+length(index1(:))-1,:)=[index1(:),index2(:),index3(:)];
    pos = pos+length(index1(:));
end
[index1,index2]=meshgrid(1:size(DofIndex_sparse,1));
Table_sparse = [index1(:),index2(:)];
Table_sparse((BaseFunIndex_x(DofIndex_sparse(Table_sparse(:,1),1),5)<=BaseFunIndex_x(DofIndex_sparse(Table_sparse(:,2),1),4))...
    |(BaseFunIndex_x(DofIndex_sparse(Table_sparse(:,1),1),4)>=BaseFunIndex_x(DofIndex_sparse(Table_sparse(:,2),1),5))...
    |(BaseFunIndex_y(DofIndex_sparse(Table_sparse(:,1),2),5)<=BaseFunIndex_y(DofIndex_sparse(Table_sparse(:,2),2),4))...
    |(BaseFunIndex_y(DofIndex_sparse(Table_sparse(:,1),2),4)>=BaseFunIndex_y(DofIndex_sparse(Table_sparse(:,2),2),5))...
    |(BaseFunIndex_z(DofIndex_sparse(Table_sparse(:,1),3),5)<=BaseFunIndex_z(DofIndex_sparse(Table_sparse(:,2),3),4))...
    |(BaseFunIndex_z(DofIndex_sparse(Table_sparse(:,1),3),4)>=BaseFunIndex_z(DofIndex_sparse(Table_sparse(:,2),3),5)),:)=[];
end

