function j = MaxJ2Sparse(jcoef,j0,J)
%MAXJ2SPARSE 根据每个维度的系数计算最大的j
% j1coef*j1+j2coef*j2+j3coef*j3<=J
% jcoef：3维向量，分别为各维度的稀疏网格参数

[~,jcoef_min]=min(jcoef);
jcoef_other=[1;2;3];
jcoef_other(jcoef_min)=[];
j=J-(jcoef(jcoef_other(1))+jcoef(jcoef_other(2)))*j0;
end

