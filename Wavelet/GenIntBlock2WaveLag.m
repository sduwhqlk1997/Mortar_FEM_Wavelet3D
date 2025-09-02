function [IntegralBlockNNZ1,IntegralBlockNNZ2] = GenIntBlock2WaveLag(NNZIndexWaveLag,DofIndex_sparse,...
    LagMeshIndex,IntegralBlock1,IntegralBlock2,DimFixed)
%GENINTBLOCK2WAVELAG 生成每个非零元的积分区间段
%   NNZIndexWaveLag：
%       1.稀疏基编号（与WaveIndex4Lag或DofIndex_sparse的行对应）
%       2.Lagrange网格编号（对应LagMeshIndex的行）
%   DofIndex_sparse:稀疏基索引，1~3列分别为对应x,y,z方向一维小波基的编号，与BaseFunIndex对应
%   LagMeshIndex：Lagrange乘子网格索引，每一行表示一个网格
%       1~2：网格第一个维度的区间端点，3~4：网格第二个维度的区间端点
%       剩下元素为对的基函数编号，即对应于LagBaseIndex的某些行,若为0则不表示任何基函数
%   IntegralBlock1,IntegralBlock2：第一和第二个维度的积分区间段
%   IntegralBlockNNZ1,IntegralBlockNNZ2：两个维度的积分区间段，每行对应一个非零元素，即与NNZIndexWaveLag的行对应

NumNNZ=size(NNZIndexWaveLag,1);
DofIndex_sparse(:,DimFixed)=[]; % 删掉固定维度
IntegralBlockNNZ1=NaN(NumNNZ,size(IntegralBlock1,2));
IntegralBlockNNZ2=NaN(NumNNZ,size(IntegralBlock2,2));
for i=1:NumNNZ
    iWave=NNZIndexWaveLag(i,1);
    iLag=NNZIndexWaveLag(i,2);
    iWave1=DofIndex_sparse(iWave,1); % 第一个维度对应的一维基函数编号，与IntegralBlock1的行对应
    iWave2=DofIndex_sparse(iWave,2); % 第二个维度对应的一维基函数编号，与IntegralBlock2的行对应
    % 第一个维度
    IntBlock=ExtractBlock(LagMeshIndex(iLag,1:2),IntegralBlock1(iWave1,:));
    IntegralBlockNNZ1(i,1:length(IntBlock))=IntBlock;
    % 第二个维度
    IntBlock=ExtractBlock(LagMeshIndex(iLag,3:4),IntegralBlock2(iWave2,:));
    IntegralBlockNNZ2(i,1:length(IntBlock))=IntBlock;
end
    function IntBlock=ExtractBlock(interval,Block) % 截取Block中落在interval区间中的区间段
        Block(isnan(Block))=[]; % 剔除缺省值
        left=Block<=interval(1);
        right=Block>=interval(2);
        ileft=any(left); % 若区间左端点不在任何区间段内则为0
        iright=any(right); % 若区间右端点不在任何区间段内则为0
        Block(left|right)=[]; % 删掉不在区间中的元素
        if  ileft&&iright  % 区间在区间段内
            IntBlock=[interval(1),Block,interval(2)];
        elseif ~ileft&&iright
            IntBlock=[Block,interval(2)];
        elseif ileft&&~iright
            IntBlock=[interval(1),Block];
        elseif ~ileft&&~iright
            IntBlock=Block;
        end
        
    end
end

