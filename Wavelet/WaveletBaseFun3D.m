function y = WaveletBaseFun3D(X1,X2,X3,interval,j,k,diff,BaseType,SpaceType)
%WAVELETBASEFUN3D 一般长方体区域上三维的小波基函数
%   X1,X2,X3:要计算的点的三个维度
%   interval：定义域，为3x2矩阵，每一行表示一个维度
%   j：三个方向小波基的level，3维向量
%   k：三个方向小波基的编号，3维向量
%   diff：三个方向的导数阶，3维向量
%   BaseType：基函数类型
%   SpaceType：函数空间类型，3维向量，0为尺度基，1为小波基

y=WaveletBaseFun1D(X1,interval(1,:),j(1),k(1),diff(1),BaseType,SpaceType(1)).*...
    WaveletBaseFun1D(X2,interval(2,:),j(2),k(2),diff(2),BaseType,SpaceType(2)).*...
    WaveletBaseFun1D(X3,interval(3,:),j(3),k(3),diff(3),BaseType,SpaceType(3));
end

