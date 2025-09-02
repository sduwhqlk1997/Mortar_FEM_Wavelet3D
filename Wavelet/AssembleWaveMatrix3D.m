function K = AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K_x,K_y,K_z)
%AssembleWaveMatrix3D 由一维刚度矩阵组装三维稀疏基刚度矩阵
%   K: 二维稀疏基对应的刚度矩阵
%   DofIndex_sparse:稀疏基索引，1~3列分别为对应x,y,z方向一维小波基的编号
%   Table_sparse 稀疏基刚度矩阵的非0元素索引
%   K_x、K_y、K_z: x、y和z方向的一维刚度矩阵
    %比如要计算diff_test=[1,0,0]; diff_trail=[0,0,1]的刚度矩阵，这时
    %K_x=AssembleMatrix1D(Coeff,interval,Table,supp,1,0,type_test,type_trail,order)
    %K_y=AssembleMatrix1D(Coeff,interval,Table,supp,0,0,type_test,type_trail,order)
    %K_z=AssembleMatrix1D(Coeff,interval,Table,supp,0,1,type_test,type_trail,order)

N = size(DofIndex_sparse,1);
ii = Table_sparse(:,1); % 行下标
jj = Table_sparse(:,2); % 列下标
linearInd_x = sub2ind(size(K_x),DofIndex_sparse(Table_sparse(:,1),1),...
    DofIndex_sparse(Table_sparse(:,2),1)); % 与x相关的积分
linearInd_y = sub2ind(size(K_y),DofIndex_sparse(Table_sparse(:,1),2),...
    DofIndex_sparse(Table_sparse(:,2),2)); % 与y相关的积分
linearInd_z = sub2ind(size(K_z),DofIndex_sparse(Table_sparse(:,1),3),...
    DofIndex_sparse(Table_sparse(:,2),3)); % 与z相关的积分
kk = K_x(linearInd_x).*K_y(linearInd_y).*K_z(linearInd_z);
K = sparse(ii,jj,kk,N,N);
end

