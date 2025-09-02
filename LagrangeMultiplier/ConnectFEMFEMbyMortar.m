function [MatMasterSider,MatSlaveSider] = ConnectFEMFEMbyMortar(MasterSiderInfo,SlaveSiderInfo,order)
%CONNECTFEMFEMBYMORTAR 使用mortar方法拼接两个使用有限元离散的区域
%   MasterSiderInfo,SlaveSiderInfo：MasterSider和SlaveSide的信息表，为cell数组：
%       1：区域，为3x2矩阵，一行表示一个维度的区间
%       2：基函数类型，为字符串
%       3：x方向的网格节点
%       4：y方向的网格节点
%       5：z方向的网格节点
%       6：有限元网格点，即Pb
%       7：有限元网格，即Tb
%   order：Gauss积分阶数，不输入默认为4阶

if ~exist('order','var') || isempty(order)
    order=4;
end
% 确定固定的维度
intervalMaster=cell2mat(MasterSiderInfo(1));
intervalSlave=cell2mat(SlaveSiderInfo(1));
[DimFixed,NodeFixed]=find([intervalMaster(:,1)==intervalSlave(:,2),...
    intervalMaster(:,2)==intervalSlave(:,1)]);
NodeFixed=intervalMaster(DimFixed,NodeFixed);
DimFree=[1;2;3];
DimFree(DimFixed)=[];
% if isempty(DimFixed)
%     exit('两区域不相接')
% end
% 生成Lagrange乘子信息矩阵
BaseTypeLag=cell2mat(SlaveSiderInfo(2));
x1StepLag=cell2mat(SlaveSiderInfo(2+DimFree(1)));
x2StepLag=cell2mat(SlaveSiderInfo(2+DimFree(2)));
[BaseIndexLag,MeshIndexLag] = GenLagBaseIndex(BaseTypeLag,x1StepLag,x2StepLag);
intervalLag([DimFixed;DimFree],:)=[NodeFixed,NodeFixed;intervalSlave(DimFree,:)];
% 提取两侧交界面与Lagrange网格相交的单元
% MasterSider
[RegionBaseIndex2DMaster,RegionMeshIndex2DMaster] ...
    = Extract2DRegionFEM(cell2mat(MasterSiderInfo(6)),cell2mat(MasterSiderInfo(7)),intervalLag);
% SlaveSider
[RegionBaseIndex2DSlave,RegionMeshIndex2DSlave] ...
    = Extract2DRegionFEM(cell2mat(SlaveSiderInfo(6)),cell2mat(SlaveSiderInfo(7)),intervalLag);
% 组装矩阵
MatMasterSider = ...
    AssembleFEMLagMat(RegionBaseIndex2DMaster,RegionMeshIndex2DMaster,...
    BaseIndexLag,MeshIndexLag,cell2mat(MasterSiderInfo(6)),cell2mat(MasterSiderInfo(7)),...
    DimFixed,NodeFixed,cell2mat(MasterSiderInfo(2)),BaseTypeLag,order);
MatSlaveSider = ...
    AssembleFEMLagMat(RegionBaseIndex2DSlave,RegionMeshIndex2DSlave,...
    BaseIndexLag,MeshIndexLag,cell2mat(SlaveSiderInfo(6)),cell2mat(SlaveSiderInfo(7)),...
    DimFixed,NodeFixed,cell2mat(SlaveSiderInfo(2)),BaseTypeLag,order);
end

