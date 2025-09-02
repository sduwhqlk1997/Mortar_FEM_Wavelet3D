function value = FEMBaseFun(X1,X2,X3,baseType,ele,baseNum,diff)
%FEMBASEFUN 一般长方体单元上的有限元基函数
%   X1,X2,X3：每一个表示一个维度的坐标
%   value：函数值
%   baseType：基函数种类
%   ele：所属单元,3x2矩阵，每一行表示一个维度
%   baseNum：基函数编号，为标量
%   diff：导数阶，为3x1向量，分别为两个维度的导数阶

X1 = (2*X1-(ele(1,1)+ele(1,2)).')./(ele(1,2)-ele(1,1)).'; % 把要计算的点变化到参考单元[-1,1]^3上
X2 = (2*X2-(ele(2,1)+ele(2,2)).')./(ele(2,2)-ele(2,1)).'; % 把要计算的点变化到参考单元[-1,1]^3上
X3 = (2*X3-(ele(3,1)+ele(3,2)).')./(ele(3,2)-ele(3,1)).'; % 把要计算的点变化到参考单元[-1,1]^3上
value=FEMBaseFunRef(baseType,X1,X2,X3,baseNum,diff(1),diff(2),diff(3));
if diff(1)==1
    value=value*2/(ele(1,2)-ele(1,1));
end
if diff(2)==1
    value=value*2/(ele(2,2)-ele(2,1));
end
if diff(3)==1
    value=value*2/(ele(3,2)-ele(3,1));
end
end

