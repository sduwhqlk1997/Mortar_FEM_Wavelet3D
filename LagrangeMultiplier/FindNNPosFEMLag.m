function NNZIndexFEMLag = FindNNPosFEMLag(FEMRegionMeshIndex2D,LagMeshIndex)
%FINDNNPOSFEMLAG 寻找Lagrange乘子与FEM矩阵相交的网格
%   FEMRegionMeshIndex2D：2D区域上的单元信息表，每一行代表一个单元
%       1~2：第一个维度的区间端点；3~4：第二个维度的区间端点；
%       5：所属单元，与Tb的行对应；剩下元素为属于该单元的面上的基函数编号（对应FEMRegionBaseIndex2D的行）
%          ，若为0则不表示任何基函数
%   LagMeshIndex：Lagrange乘子网格索引，每一行表示一个网格
%       1~2：网格第一个维度的区间端点，3~4：网格第二个维度的区间端点
%       剩下元素为对的基函数编号，即对应于LagBaseIndex的某些行,若为0则不表示任何基函数
%   NNZIndexFEMLag：FEM与Lagrange乘子相交的网格，每一行表示一组相交网格
%       1：有限元网格编号（对应FEMRegionMeshIndex2D的行）
%       2：Lagrange网格编号（对应LagMeshIndex的行）
%       3~4：第一个维度的区间交集（可用于后续数值积分）
%       5~6：第二个维度的区间交集（可用于后续数值积分）

NFEMMesh=size(FEMRegionMeshIndex2D,1);
NLagMesh=size(LagMeshIndex,1);
IndexFEM=(1:NFEMMesh).';
IndexLag=(1:NLagMesh).';
[IndexFEM,IndexLag]=meshgrid(IndexFEM,IndexLag);
IndexFEM=IndexFEM(:);
IndexLag=IndexLag(:);
IndexDel=FEMRegionMeshIndex2D(IndexFEM,1)>=LagMeshIndex(IndexLag,2)|...
    FEMRegionMeshIndex2D(IndexFEM,2)<=LagMeshIndex(IndexLag,1)|...
    FEMRegionMeshIndex2D(IndexFEM,3)>=LagMeshIndex(IndexLag,4)|...
    FEMRegionMeshIndex2D(IndexFEM,4)<=LagMeshIndex(IndexLag,3); % 所有不相交的网格
IndexFEM(IndexDel)=[];
IndexLag(IndexDel)=[];
NNZIndexFEMLag=...
    [IndexFEM,IndexLag,...
    max([FEMRegionMeshIndex2D(IndexFEM,1),LagMeshIndex(IndexLag,1)],[],2),...
    min([FEMRegionMeshIndex2D(IndexFEM,2),LagMeshIndex(IndexLag,2)],[],2),...
    max([FEMRegionMeshIndex2D(IndexFEM,3),LagMeshIndex(IndexLag,3)],[],2),...
    min([FEMRegionMeshIndex2D(IndexFEM,4),LagMeshIndex(IndexLag,4)],[],2)];
end

