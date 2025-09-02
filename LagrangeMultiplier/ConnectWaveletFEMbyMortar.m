function [MatWaveSide,MatFEMSide] = ConnectWaveletFEMbyMortar(WaveSideInfo,FEMSideInfo,order)
%CONNECTWAVELETFEMBYMORTAR 连接小波和有限元离散的区域，有限元作为SlaveSide
%   WaveSideInfo：小波离散区域的信息表
%       1：区域，为3x2矩阵，一行表示一个维度的区间
%       2：基函数类型，为字符串
%       3：BaseFunIndex_x （变量说明见GenWaveBaseIndex4Lag.m）
%       4：BaseFunIndex_y
%       5：BaseFunIndex_z
%       6：DofIndex_sparse（稀疏基索引，1~3列分别为对应x,y,z方向一维小波基的编号，与BaseFunIndex对应）
%       7：IntegralBlock_x（每个一维小波基的积分区间段，详细说明见GenWaveletBaseIndex1D.m）
%       8：IntegralBlock_y
%       9：IntegralBlock_z
%   FEMSideInfo：有限元离散区域的信息表
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
intervalWave=cell2mat(WaveSideInfo(1));
intervalFEM=cell2mat(FEMSideInfo(1));
[DimFixed,NodeFixed]=find([intervalWave(:,1)==intervalFEM(:,2),...
    intervalWave(:,2)==intervalFEM(:,1)]);
NodeFixed=intervalWave(DimFixed,NodeFixed);
DimFree=[1;2;3];
DimFree(DimFixed)=[];
% 生成Lagrange乘子信息矩阵
BaseTypeLag=cell2mat(FEMSideInfo(2));
x1StepLag=cell2mat(FEMSideInfo(2+DimFree(1)));
x2StepLag=cell2mat(FEMSideInfo(2+DimFree(2)));
[BaseIndexLag,MeshIndexLag] = GenLagBaseIndex(BaseTypeLag,x1StepLag,x2StepLag);
intervalLag([DimFixed;DimFree],:)=[NodeFixed,NodeFixed;intervalFEM(DimFree,:)];
% 提取两侧交界面与Lagrange网格相交的单元
% WaveletSide
WaveIndex4Lag = GenWaveBaseIndex4Lag(cell2mat(WaveSideInfo(6)),...
    cell2mat(WaveSideInfo(3)),cell2mat(WaveSideInfo(4)),cell2mat(WaveSideInfo(5)));
WaveletRegionBaseIndex2D = Extract2DRegionWavelet(WaveIndex4Lag,intervalLag);
% FEMSide
[RegionBaseIndex2DFEM,RegionMeshIndex2DFEM] ...
    = Extract2DRegionFEM(cell2mat(FEMSideInfo(6)),cell2mat(FEMSideInfo(7)),intervalLag);
% 组装矩阵
% WaveletSide
MatWaveSide = AssembleWaveLagMat(cell2mat(WaveSideInfo(6)),MeshIndexLag,BaseIndexLag,...
    cell2mat(WaveSideInfo(6+DimFree(1))),cell2mat(WaveSideInfo(6+DimFree(2))),...
    cell2mat(WaveSideInfo(3)),cell2mat(WaveSideInfo(4)),cell2mat(WaveSideInfo(5)),...
    WaveletRegionBaseIndex2D,WaveIndex4Lag,...
    DimFixed,NodeFixed,cell2mat(WaveSideInfo(2)),BaseTypeLag,order,intervalWave);
% FEMSide
MatFEMSide = ...
    AssembleFEMLagMat(RegionBaseIndex2DFEM,RegionMeshIndex2DFEM,...
    BaseIndexLag,MeshIndexLag,cell2mat(FEMSideInfo(6)),cell2mat(FEMSideInfo(7)),...
    DimFixed,NodeFixed,cell2mat(FEMSideInfo(2)),BaseTypeLag,order);
end

