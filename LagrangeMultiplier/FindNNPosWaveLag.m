function NNZIndexWaveLag = FindNNPosWaveLag(WaveletRegionBaseIndex2D,...
    LagMeshIndex,WaveIndex4Lag,DimFixed)
%FINDNNPOSWAVELAG 此处显示有关此函数的摘要
%   WaveletRegionBaseIndex2D：对应的稀疏基编号，与DofIndex_sparse的行对应(或者说WaveIndex4Lag的行)
%   LagMeshIndex：Lagrange乘子网格索引，每一行表示一个网格
%       1~2：网格第一个维度的区间端点，3~4：网格第二个维度的区间端点
%       剩下元素为对的基函数编号，即对应于LagBaseIndex的某些行,若为0则不表示任何基函数
%   WaveIndex4Lag：用于组装Lagrange-Wavelet连接矩阵的索引矩阵(每一行与DofIndex_sparse的行对应)
%       1-2：x方向的支集端点；3-4：y方向的支集端点；5-6：z方向的支集端点
%       7-9：对应三个方向基函数的编号，分别与BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z的行对应
%   DimFixed：固定维度

%   NNZIndexWaveLag：
%       1.稀疏基编号（与WaveIndex4Lag或DofIndex_sparse的行对应）
%       2.Lagrange网格编号（对应LagMeshIndex的行）
%   IntegralBlockNNZ1,IntegralBlockNNZ2：两个维度的积分区间段，每行对应一个非零元素，即与NNZIndexWaveLag的行对应
Nwave=length(WaveletRegionBaseIndex2D);
NLagMesh=size(LagMeshIndex,1);
IndexWave=(1:Nwave).';
IndexLag=(1:NLagMesh).';
[IndexWave,IndexLag]=meshgrid(IndexWave,IndexLag);
IndexWave=IndexWave(:);
IndexLag=IndexLag(:);
WaveIndex4Lag=WaveIndex4Lag(WaveletRegionBaseIndex2D,:); % 提取在交界区域的稀疏基信息
WaveIndex4Lag(:,[2*DimFixed-1,2*DimFixed,6+DimFixed])=[]; 
%{此时WaveIndex4Lag各列含义变为
%1-2：第一个维度支集端点；3-4：第二个维度支集端点
%5-6：对应两个维度的基函数编号，
%}
IndexDel=WaveIndex4Lag(IndexWave,1)>=LagMeshIndex(IndexLag,2)|...
    WaveIndex4Lag(IndexWave,2)<=LagMeshIndex(IndexLag,1)|...
    WaveIndex4Lag(IndexWave,3)>=LagMeshIndex(IndexLag,4)|...
    WaveIndex4Lag(IndexWave,4)<=LagMeshIndex(IndexLag,3);
IndexWave(IndexDel)=[];
IndexLag(IndexDel)=[];
NNZIndexWaveLag=[WaveletRegionBaseIndex2D(IndexWave),IndexLag];
end

