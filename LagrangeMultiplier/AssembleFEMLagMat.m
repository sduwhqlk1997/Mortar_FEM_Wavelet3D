function K = ...
    AssembleFEMLagMat(FEMRegionBaseIndex2D,FEMRegionMeshIndex2D,...
    LagBaseIndex,LagMeshIndex,Pb,Tb,DimFixed,NodeFixed,baseTypeFEM,baseTypeLag,order)
%{ASSEMBLEFEMLAGMAT 组装Lagrange乘子矩阵块（右上角块，即：行为有限元坐标，列为lagrange坐标）
   % FEMRegionBaseIndex2D：2D区域上的有限元基函数信息表，每一行代表一个基函数
   %     1~2：第一个维度的区间端点；3~4：第二个维度的区间端点；
   %     5：所属单元，与Tb的行对应；6：节点相对编号（对应单元中的序号）7：节点绝对编号（对应Pb中的编号）
   % FEMRegionMeshIndex2D：2D区域上的有限元单元信息表，每一行代表一个单元
   %     1~2：第一个维度的区间端点；3~4：第二个维度的区间端点；
   %     5：所属单元，与Tb的行对应；剩下元素为属于该单元的面上的基函数编号，若为0则不表示任何基函数
   % LagBaseIndex：Lagrange乘子基函数信息表，每一行表示一个基函数，各列存储的信息如下所示：
   %     1~4：所属单元，1~2第一个维度，3~4第二个维度
   %     5~6：两个维度对应的一维基函数编号
   %     7~8：标记两个维度是否为边界单元，若为0则为内部单元，若为1则为左侧单元，若为2则为右侧单元
   % LagMeshIndex：Lagrange乘子网格索引，每一行表示一个网格
   %     1~2：网格第一个维度的区间，3~4：网格第二个维度的区间
   %     剩下元素为对的基函数编号，即对应于LagBaseIndex的某些行,若为0则不表示任何基函数
   % Pb：总区域的有限元节点
   % Tb：总区域的有限元网格
   % DimFixed：固定的维度，1：x方向，2：y方向，3：z方向
   % NodeFixed：固定维度的固定点，为标量
   % baseType：基函数类型
   % order：Gauss积分的阶数
   %}

NFEMDOF=size(Pb,1); % 有限元自由度总数，即为K的行数
NNZIndexFEMLag = FindNNPosFEMLag(FEMRegionMeshIndex2D,LagMeshIndex);
Triplets=[];
for i=1:size(NNZIndexFEMLag,1)
    BaseListFEM=FEMRegionMeshIndex2D(NNZIndexFEMLag(i,1),6:end);
    BaseListLag=LagMeshIndex(NNZIndexFEMLag(i,2),5:end);
    BaseListFEM(BaseListFEM==0)=[];
    BaseListLag(BaseListLag==0)=[];
    % 生成单元所属的区间
    NumFEMele=FEMRegionMeshIndex2D(NNZIndexFEMLag(i,1),5); % FEM单元编号，对应Tb的行
    FEMele=[min(Pb(Tb(NumFEMele,:),1)),max(Pb(Tb(NumFEMele,:),1));...
    min(Pb(Tb(NumFEMele,:),2)),max(Pb(Tb(NumFEMele,:),2));...
    min(Pb(Tb(NumFEMele,:),3)),max(Pb(Tb(NumFEMele,:),3))];
    Lagele=[LagMeshIndex(NNZIndexFEMLag(i,2),1:2);...
        LagMeshIndex(NNZIndexFEMLag(i,2),3:4)];
    
    for iFEM=BaseListFEM % 对应FEMRegionBaseIndex2D的行
        for iLag=BaseListLag % 对应LagBaseIndex的行
            switch DimFixed % 要计算积分的有限元基函数
                case 1
                    FunFEM=@(X,Y) FEMBaseFun(NodeFixed*ones(size(X,1),size(X,2)),X,Y,...
                        baseTypeFEM,...
                        FEMele,...
                        FEMRegionBaseIndex2D(iFEM,6),[0;0;0]);
                case 2
                    FunFEM=@(X,Y) FEMBaseFun(X,NodeFixed*ones(size(X,1),size(X,2)),Y,...
                        baseTypeFEM,...
                        FEMele,...
                        FEMRegionBaseIndex2D(iFEM,6),[0;0;0]);
                case 3
                    FunFEM=@(X,Y) FEMBaseFun(X,Y,NodeFixed*ones(size(X,1),size(X,2)),...
                        baseTypeFEM,...
                        FEMele,...
                        FEMRegionBaseIndex2D(iFEM,6),[0;0;0]);
            end
            FunLag=@(X,Y) LagBaseFun(X,Y,LagBaseIndex(iLag,7:8)',...
                baseTypeLag,Lagele,LagBaseIndex(iLag,5:6)',[0;0]);
            f=@(X,Y)FunFEM(X,Y).*FunLag(X,Y);
            value=GaussIntegralBlocks(f,order,2,NNZIndexFEMLag(i,3:4)',NNZIndexFEMLag(i,5:6)');
            Triplets=[Triplets;FEMRegionBaseIndex2D(iFEM,7),iLag,value];
        end
    end
end
K=sparse(Triplets(:,1),Triplets(:,2),Triplets(:,3),NFEMDOF,size(LagBaseIndex,1));
end

