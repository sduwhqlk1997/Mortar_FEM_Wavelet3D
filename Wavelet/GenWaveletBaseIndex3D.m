function [BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z,...
    IntegralBlock_x,IntegralBlock_y,IntegralBlock_z] ...
    = GenWaveletBaseIndex3D(interval,j0,j,BaseType)
%GENWAVELETBASEINDEX3D 此处显示有关此函数的摘要
%   BaseFunIndex：一维小波基函数信息表
%       1：SpaceType，0表示尺度基，1表示小波基
%       2：level
%       3：基函数编号
%       4：支集左端点
%       5：支集右端点
%       6：基函数坐标索引
%       7：基函数编号
%   IntegralBlock：一维小波基函数的积分分段，为矩阵，每一行表示一个基函数，非NaN数据表示节点（从左到右依次排列）
%       顺序与BaseFunIndex一致

[BaseFunIndex_ref,IntegralBlock_ref] = GenWaveletBaseIndex1D([0,1],j0,j,BaseType);
BaseFunIndex_x=BaseFunIndex_ref;
IntegralBlock_x=IntegralBlock_ref;
BaseFunIndex_x(:,4:6) = (interval(1,2)-interval(1,1)).*BaseFunIndex_x(:,4:6)+interval(1,1);
IntegralBlock_x=(interval(1,2)-interval(1,1)).*IntegralBlock_x+interval(1,1);

BaseFunIndex_y=BaseFunIndex_ref;
IntegralBlock_y=IntegralBlock_ref;
BaseFunIndex_y(:,4:6) = (interval(2,2)-interval(2,1)).*BaseFunIndex_y(:,4:6)+interval(2,1);
IntegralBlock_y=(interval(2,2)-interval(2,1)).*IntegralBlock_y+interval(2,1);

BaseFunIndex_z=BaseFunIndex_ref;
IntegralBlock_z=IntegralBlock_ref;
BaseFunIndex_z(:,4:6) = (interval(3,2)-interval(3,1)).*BaseFunIndex_z(:,4:6)+interval(3,1);
IntegralBlock_z=(interval(3,2)-interval(3,1)).*IntegralBlock_z+interval(3,1);
end

