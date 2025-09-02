function [X,Y,Z] = GenWaveGaussBlock(interval_Sol,interval,dim,BaseIndex_x,DofX_num,BaseIndex_y,DofY_num,BaseIndex_z,DofZ_num)
%GenWaveGaussBlock 生成分片Gauss积分的分片信息
%   interval_Sol:问题的求解域
%   interval:积分区间，为dim*2的矩阵，每一行代表一个维度
%   dim:积分维度
%   BaseIndex_x,BaseIndex_y,BaseIndex_z:三个坐标方向的基函数索引，每一行代表一个基函数，其中第二列表示该基函数的level
%   DofX_num,DofY_num,DofZ_num:各方向基函数的自由度编号为列向量
%   X,Y,Z:三个方向的积分分片
if dim==1
    jmax = max(BaseIndex_x(DofX_num,2)); % 获取该维度的最高level
    d = 2^(-jmax-1)*(interval_Sol(1,2)-interval_Sol(1,1));
    X = (interval(1,1):d:interval(1,2))';
    if length(X)<=1
        X = interval;
    end
elseif dim==2
    jmax_x = max(BaseIndex_x(DofX_num,2));
    jmax_y = max(BaseIndex_y(DofY_num,2));
    
    d = 2^(-jmax_x-1)*(interval_Sol(1,2)-interval_Sol(1,1));
    X = (interval(1,1):d:interval(1,2))';
    if length(X)<=1
        X = interval(1,:);
    end
    
    d = 2^(-jmax_y-1)*(interval_Sol(2,2)-interval_Sol(2,1));
    Y = (interval(2,1):d:interval(2,2))';
    if length(Y)<=1
        Y = interval(2,:);
    end
elseif dim==3 
    jmax_x = max(BaseIndex_x(DofX_num,2));
    jmax_y = max(BaseIndex_y(DofY_num,2));
    jmax_z = max(BaseIndex_z(DofZ_num,2));
    
    d = 2^(-jmax_x-1)*(interval_Sol(1,2)-interval_Sol(1,1));
    X = (interval(1,1):d:interval(1,2))';
    if length(X)<=1
        % 1
        d = 2^(-jmax_x-2)*(interval_Sol(1,2)-interval_Sol(1,1));
        X = (interval(1,1):d:interval(1,2))';
    end
    
    d = 2^(-jmax_y-1)*(interval_Sol(2,2)-interval_Sol(2,1));
    Y = (interval(2,1):d:interval(2,2))';
    if length(Y)<=1
        % 1
        d = 2^(-jmax_y-2)*(interval_Sol(2,2)-interval_Sol(2,1));
        Y = (interval(2,1):d:interval(2,2))';
    end
    
    d = 2^(-jmax_z-1)*(interval_Sol(3,2)-interval_Sol(3,1));
    Z = (interval(3,1):d:interval(3,2))';
    if length(Z)<=1
        % 1
        d = 2^(-jmax_z-2)*(interval_Sol(3,2)-interval_Sol(3,1));
        Z = (interval(3,1):d:interval(3,2))';
    end
else
    error("维度不能大于三维，不能小于一维")
end
end

