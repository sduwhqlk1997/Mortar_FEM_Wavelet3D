function y = WaveletBaseFun1D(x,interval,j,k,diff,BaseType,SpaceType)
%WAVELETBASEFUN1D 一般区间interval上的小波基
%   x:自变量的值
%   interval:定义域
%   j:基函数的level
%   k:基函数编号
%   diff: 导数阶
%   BaseType: 基函数类型
%   SpaceType：函数空间类型，0为尺度基，1为小波基

x = (x-interval(1))./(interval(2)-interval(1)); % 将要计算的点映射到参考单元上
switch diff
    case 0
        y = WaveletBaseRef(BaseType,SpaceType,x,j,k,diff);
    case 1
        y = WaveletBaseRef(BaseType,SpaceType,x,j,k,diff)./(interval(2)-interval(1));
    otherwise
        error("导数阶数只能为0阶或1阶")
end
end

