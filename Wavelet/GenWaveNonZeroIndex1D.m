function Table = GenWaveNonZeroIndex1D(supp)
%GenWaveNonZeroIndex1D 生成紧支集交集非空的一维基函数索引组合
%   Table: 非零元素索引表 1~4个元素依次为：第一个基函数编号、第二个基函数编号、紧支集交集左端点、紧支集交集右端点
%   supp: 每行7个元素，依次为type2(尺度基(0)or小波基(1))、j、k、支集左端点、支集右端点、基函数的坐标索引和基函数编号 
[index1,index2]=meshgrid(supp(:,end));
index = [index1(:),index2(:)];
index(index(:,1)>index(:,2),:)=[];
index((supp(index(:,1),5)<=supp(index(:,2),4))|(supp(index(:,1),4)>=supp(index(:,2),5)),:)=[];
Table = [index,max(supp(index(:,1),4),supp(index(:,2),4)),min(supp(index(:,1),5),supp(index(:,2),5))];
end

