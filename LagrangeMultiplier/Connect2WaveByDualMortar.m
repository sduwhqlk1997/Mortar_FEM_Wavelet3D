function [MatMaster,MatSlave] = Connect2WaveByDualMortar(MasterSideInfo,SlaveSideInfo)
%CONNECT2WAVEBYDUALMORTAR 用小波基的对偶基来连接两个区域，要求：
% 1.两个区域用同一种小波离散，可以不同level
% 2.两个区域的端点必须重合
%   DimFixed：固定的维度
%   NodeFixed：固定维度的坐标
%   MasterSideInfo,SlaveSideInfo：两侧区域小波基的信息表，为cell数组，各元素含义以此为：
%       1.区域，3x2矩阵，每一行表示一个维度
%       2.三个方向的level（j）
%       3.三个方向的index（k）
%       4.三个方向的SpaceType，0：尺度基，1：小波基
%       5.三个方向的坐标索引
%       6.自由度类型

% 确定固定维度
intervalMaster=cell2mat(MasterSideInfo(1));
intervalSlave=cell2mat(SlaveSideInfo(1));
[DimFixed,NodeFixed]=find([intervalMaster(:,1)==intervalSlave(:,2),...
    intervalMaster(:,2)==intervalSlave(:,1)]);
NodeFixed=intervalMaster(DimFixed,NodeFixed);
DimFree=[1;2;3];
DimFree(DimFixed)=[];
% 提取交界面处的自由度
CoorIndexMaster=cell2mat(MasterSideInfo(5));
CoorIndexSlave=cell2mat(SlaveSideInfo(5));
InterfaceIndexMaster=CoorIndexMaster(:,DimFixed)==NodeFixed;
InterfaceIndexSlave=CoorIndexSlave(:,DimFixed)==NodeFixed;
% 构造Lagrange乘子基信息表
% LagIndex=InterfaceIndexSlave&...
%     CoorIndexSlave(:,DimFree(1))~=intervalSlave(DimFree(1),1)&...
%     CoorIndexSlave(:,DimFree(1))~=intervalSlave(DimFree(1),2)&...
%     CoorIndexSlave(:,DimFree(2))~=intervalSlave(DimFree(2),1)&...
%     CoorIndexSlave(:,DimFree(2))~=intervalSlave(DimFree(2),2); % mortar元不需要边界自由度
LagIndex=InterfaceIndexSlave;
jSlave=cell2mat(SlaveSideInfo(2));
kSlave=cell2mat(SlaveSideInfo(3));
SpaceTypeSlave=cell2mat(SlaveSideInfo(4));
DofTypeSlave=cell2mat(SlaveSideInfo(6));
jLag=jSlave(LagIndex,DimFree);
kLag=kSlave(LagIndex,DimFree);
SpaceTypeLag=SpaceTypeSlave(LagIndex,DimFree);
DofTypeLag = DofTypeSlave(LagIndex);
% 组装MasterSide的连接矩阵
jMaster=cell2mat(MasterSideInfo(2));
kMaster=cell2mat(MasterSideInfo(3));
SpaceTypeMaster=cell2mat(MasterSideInfo(4));
DofTypeMaster=cell2mat(MasterSideInfo(6));
InterfaceIndexMaster=find(InterfaceIndexMaster);
[MasterIndex,LagIndex]=find(SpaceTypeMaster(InterfaceIndexMaster,DimFree(1))==SpaceTypeLag(:,1)'&...
    SpaceTypeMaster(InterfaceIndexMaster,DimFree(2))==SpaceTypeLag(:,2)'&...
    jMaster(InterfaceIndexMaster,DimFree(1))==jLag(:,1)'&...
    jMaster(InterfaceIndexMaster,DimFree(2))==jLag(:,2)'&...
    kMaster(InterfaceIndexMaster,DimFree(1))==kLag(:,1)'&...
    kMaster(InterfaceIndexMaster,DimFree(2))==kLag(:,2)'&...
    DofTypeMaster(InterfaceIndexMaster)==DofTypeLag');
MasterIndex=InterfaceIndexMaster(MasterIndex);
MatMaster=sparse(MasterIndex,LagIndex,ones(length(MasterIndex),1),size(jMaster,1),size(jLag,1));
% 组装SlaveSide的连接矩阵
InterfaceIndexSlave=find(InterfaceIndexSlave);
[SlaveIndex,LagIndex]=find(SpaceTypeSlave(InterfaceIndexSlave,DimFree(1))==SpaceTypeLag(:,1)'&...
    SpaceTypeSlave(InterfaceIndexSlave,DimFree(2))==SpaceTypeLag(:,2)'&...
    jSlave(InterfaceIndexSlave,DimFree(1))==jLag(:,1)'&...
    jSlave(InterfaceIndexSlave,DimFree(2))==jLag(:,2)'&...
    kSlave(InterfaceIndexSlave,DimFree(1))==kLag(:,1)'&...
    kSlave(InterfaceIndexSlave,DimFree(2))==kLag(:,2)'&...
    DofTypeSlave(InterfaceIndexSlave)==DofTypeLag');
SlaveIndex=InterfaceIndexSlave(SlaveIndex);
MatSlave=sparse(SlaveIndex,LagIndex,ones(length(SlaveIndex),1),size(jSlave,1),size(jLag,1));
end

