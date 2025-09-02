function [FEMRegionBaseIndex2D,FEMRegionMeshIndex2D] = Extract2DRegionFEM(Pb,Tb,interval)
%EXTRACT2DREGIONFEM 从有限元网格和节点信息中提取两维区域的信息，并生成信息表
%   baseType：基函数类型
%   Pb,Tb：网格点和网格信息
%   interval：3x2矩阵，每一行表示一个维度
%   FEMRegionBaseIndex2D：2D区域上的基函数信息表，每一行代表一个基函数
%       1~2：第一个维度的区间端点；3~4：第二个维度的区间端点；
%       5：所属单元，与Tb的行对应；6：节点相对编号（对应单元中的序号）7：节点绝对编号（对应Pb中的编号）
%   FEMRegionMeshIndex2D：2D区域上的单元信息表，每一行代表一个单元
%       1~2：第一个维度的区间端点；3~4：第二个维度的区间端点；
%       5：所属单元，与Tb的行对应；剩下元素为属于该单元的面上的基函数编号，若为0则不表示任何基函数

% 提取出interval中的单元
dimFree=[1;2;3];
dimFixed=find(interval(:,1)==interval(:,2)); % 固定的维度
dimFree(dimFixed)=[];
dimFree2TbInterval=2*dimFree-1; % 每个维度起始点对应TbInterval中的位置
dimFixed2TbInterval=2*dimFixed-1;
NumNodeElem=size(Tb,2); % 每个单元中结点数量
TbInterval=[min(reshape(Pb(Tb,1),[],NumNodeElem),[],2),max(reshape(Pb(Tb,1),[],NumNodeElem),[],2),...
    min(reshape(Pb(Tb,2),[],NumNodeElem),[],2),max(reshape(Pb(Tb,2),[],NumNodeElem),[],2),...
    min(reshape(Pb(Tb,3),[],NumNodeElem),[],2),max(reshape(Pb(Tb,3),[],NumNodeElem),[],2)];
ElemIndex=find(~(any(TbInterval(:,dimFree2TbInterval)>=interval(dimFree,2)',2)...
    | any(TbInterval(:,dimFree2TbInterval+1)<=interval(dimFree,1)',2))&...
    any(TbInterval(:,[dimFixed2TbInterval,dimFixed2TbInterval+1])==interval(dimFixed,:),2)); % 区域内所有的单元的编号
BaseIndex_ele=abs(reshape(Pb(Tb(ElemIndex,:),dimFixed),[],...
    NumNodeElem)-interval(dimFixed,1))/(interval(dimFree(1),2)-interval(dimFree(1),1))<1e-10;
%BaseIndex_ele=ElemIndex(BaseIndex_ele);
Base_pos=1; % 指向RegionBaseIndex2D的尾部
FEMRegionBaseIndex2D=zeros(nnz(BaseIndex_ele),7);
FEMRegionMeshIndex2D=zeros(length(ElemIndex),5+NumNodeElem);
for i=1:size(BaseIndex_ele,1)
    index_node=find(BaseIndex_ele(i,:));
    Nnode=length(index_node);
    FEMRegionBaseIndex2D(Base_pos:Base_pos+Nnode-1,:)=[repmat(TbInterval(ElemIndex(i),...
        [dimFree2TbInterval(1),dimFree2TbInterval(1)+1,dimFree2TbInterval(2),dimFree2TbInterval(2)+1]),Nnode,1),...
        repmat(ElemIndex(i),Nnode,1),index_node',Tb(ElemIndex(i),index_node)'];
    FEMRegionMeshIndex2D(i,1:5+Nnode)=[TbInterval(ElemIndex(i),...
        [dimFree2TbInterval(1),dimFree2TbInterval(1)+1,dimFree2TbInterval(2),dimFree2TbInterval(2)+1]),...
        ElemIndex(i),Base_pos:Base_pos+Nnode-1];
    Base_pos=Base_pos+Nnode;
end
% ElemIndex2Base=repmat(ElemIndex',NumNodeElem,1);
% ElemIndex2Base=ElemIndex2Base(:);
% %RegionBaseIndex2D=zeros(NumNodeElem*length(ElemIndex),7);
% RegionBaseIndex2D=[TbInterval(ElemIndex2Base,[dimFree2TbInterval(1),dimFree2TbInterval(1)+1,...
%     dimFree2TbInterval(2),dimFree2TbInterval(2)+1]),ElemIndex2Base,...
%     repmat((1:NumNodeElem)',length(ElemIndex),1),reshape(Tb(ElemIndex,:)',[],1)];
% RegionMeshIndex2D=[TbInterval(ElemIndex,[dimFree2TbInterval(1),dimFree2TbInterval(1)+1,...
%     dimFree2TbInterval(2),dimFree2TbInterval(2)+1]),ElemIndex,...
%     reshape(1:size(RegionBaseIndex2D,1),NumNodeElem,[]).'];
end

