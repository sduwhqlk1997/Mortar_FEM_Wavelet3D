function WaveIndex4Lag = GenWaveBaseIndex4Lag(DofIndex_sparse,...
    BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z)
%GENWAVEBASEINDEX4LAG 生成用于组装Lagrange-Wavelet矩阵的索引
%   DofIndex_sparse:稀疏基索引，1~3列分别为对应x,y,z方向一维小波基的编号，与BaseFunIndex对应
%   BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z:分别为x,y,z方向的一维参考基的基函数信息，
%       1：SpaceType，0表示尺度基，1表示小波基
%       2：level
%       3：基函数编号
%       4：支集左端点
%       5：支集右端点
%       6：基函数坐标索引
%       7：基函数编号
%   WaveIndex4Lag：用于组装Lagrange-Wavelet连接矩阵的索引矩阵(每一行与DofIndex_sparse的行对应)
%       1-2：x方向的支集端点；3-4：y方向的支集端点；5-6：z方向的支集端点
%       7-9：对应三个方向基函数的编号，分别与BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z的行对应
%       10-12：对应三个方向基函数的坐标索引

WaveIndex4Lag=[BaseFunIndex_x(DofIndex_sparse(:,1),4:5),...
    BaseFunIndex_y(DofIndex_sparse(:,2),4:5),...
    BaseFunIndex_z(DofIndex_sparse(:,3),4:5),...
    DofIndex_sparse,...
    BaseFunIndex_x(DofIndex_sparse(:,1),6),...
    BaseFunIndex_y(DofIndex_sparse(:,2),6),...
    BaseFunIndex_z(DofIndex_sparse(:,3),6)];

end

