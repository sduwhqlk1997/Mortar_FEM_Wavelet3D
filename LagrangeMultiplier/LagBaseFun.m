function value = LagBaseFun(X1,X2,flag,baseType,ele,baseNum,diff)
%LAGBASEFUN 拉格朗日乘子的基函数
%   X1,X2：要计算的点，分别表示维度1、2的坐标
%   flag：2x1向量，标识两个维度的基函数是否为最左或最右侧单元，若为0则为内部单元，
%       若为1则为左侧单元，若为2则为右侧单元
%   value：函数值
%   baseType：基函数种类
%   ele：所属单元,2x2矩阵，每一行表示一个维度
%   baseNum：二维向量，表示两个坐标方向基函数的编号
%   diff：导数阶，为2x1向量，分别为两个维度的导数阶



X1 = (2*X1-(ele(1,1)+ele(1,2)).')./(ele(1,2)-ele(1,1)).'; % 把要计算的点变化到二维参考[-1,1]x[-1,1]单元上
X2 = (2*X2-(ele(2,1)+ele(2,2)).')./(ele(2,2)-ele(2,1)).'; % 把要计算的点变化到二维参考[-1,1]x[-1,1]单元上
% X(:,1)=X(:,1)*2/(ele(1,2)-ele(1,1))-
value = LagBaseFunRef1D(X1,baseType,baseNum(1),diff(1),flag(1)).*...
    LagBaseFunRef1D(X2,baseType,baseNum(2),diff(2),flag(2));
if diff(1)==1
    value=value*2/(ele(1,2)-ele(1,1));
end
if diff(2)==1
    value=value*2/(ele(2,2)-ele(2,1));
end
end

