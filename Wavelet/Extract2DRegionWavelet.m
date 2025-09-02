function WaveletRegionBaseIndex2D = Extract2DRegionWavelet(WaveIndex4Lag,interval)
%EXTRACT2DREGIONWAVELET 在小波信息表中提取给定二维区域中基函数的信息表
%   WaveIndex4Lag：用于组装Lagrange-Wavelet连接矩阵的索引矩阵(每一行与DofIndex_sparse的行对应)
%       1-2：x方向的支集端点；3-4：y方向的支集端点；5-6：z方向的支集端点
%       7-9：对应三个方向基函数的编号，分别与BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z的行对应
%       10-12：对应三个方向基函数的坐标索引
%   interval：给定的区域
%   WaveletRegionBaseIndex2D：对应的稀疏基编号，与DofIndex_sparse的行对应(或者说WaveIndex4Lag的行)

dimFree=[1;2;3];
dimFixed=find(interval(:,1)==interval(:,2)); % 固定的维度
dimFree(dimFixed)=[];
dimFreePosStart=2*dimFree-1; % 自由的维度的区间起点对应的WaveletRegionBaseIndex2D的列
dimFreePosEnd=2*dimFree; % 自由的维度的区间终点对应的WaveletRegionBaseIndex2D的列
% 确定区域内的基函数编号
WaveletRegionBaseIndex2D=find(WaveIndex4Lag(:,9+dimFixed)==interval(dimFixed,1)&...
    ~any(WaveIndex4Lag(:,dimFreePosStart)>=interval(dimFree,2)'|...
    WaveIndex4Lag(:,dimFreePosEnd)<=interval(dimFree,1)',2));
% N=length(InterFaceIndex);
% WaveletRegionBaseIndex2D=[repmat(intervalWavelet(1,:),N,1),...
%     repmat(intervalWavelet(2,:),N,1),...
%     repmat(intervalWavelet(3,:),N,1),...
%     InterFaceIndex];
end

