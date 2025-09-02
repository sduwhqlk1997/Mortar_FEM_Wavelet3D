function value = LagBaseFunRef1D(X,baseType,baseNum,diff,flag)
%LAGBASEFUNREF1D 定义在一维参考单元[-1,1]上的Lagrange乘子基函数
%   X：一维列向量，要计算的点
%   baseType：基函数类型
%   baseNum：基函数编号
%   diff：导数阶
%   flag：是否为边界单元, 若为0则为内部单元，若为1则为左侧单元，若为2则为右侧单元
%       注：若要计算的基函数为左、右边界单元，则一定要保证该基函数不为左、右端点的基函数

switch baseType
    case "linear" % 线性基函数
        switch diff
            case 0
                f1 = @(x) (1-x)/2;
                f2 = @(x) (1+x)/2;
                switch baseNum
                    case 1
                        f=f1;
                    case 2
                        f=f2;
                end
            case 1
                f1 = @(x) -ones(length(x),1)/2;
                f2 = @(x) ones(length(x),1)/2;
                switch baseNum
                    case 1
                        f=@(x) f1;
                    case 2
                        f=@(x) f2;
                end
        end
        if flag~=0 % 如果为边界单元上的基函数
            switch flag
                case 1 % 左边界单元，此时baseNum只可能为2
                    f=@(x) f(x) + f1(x);
                case 2 % 右边界单元，此时baseNum只可能为1
                    f=@(x) f(x) + f2(x);
            end
        end
    case "quadratic" % 二次基函数
        switch diff
            case 0
                f1=@(x)x.*(x-1)/2; % 左（-1）
                f2=@(x)x.*(x+1)/2; % 右（1）
                f3=@(x)1-x.^2; % 中（0）
            case 1
                f1=@(x)x-0.5; % 左（-1）
                f2=@(x)x+0.5; % 右（1）
                f3=@(x)-2*x; % 中（0）
        end
        switch baseNum
            case 1
                f=f1;
            case 2
                f=f2;
            case 3
                f=f3;
        end
        if flag~=0 % 如果为端点上的基函数
            switch flag
                case 1  % 左端点，此时baseNum只可能为3
                    f=@(x) f(x)+f1(x);
                case 2  % 右边界单元，此时baseNum只可能为3
                    f=@(x) f(x)+f2(x);
            end
        end
end
value = f(X);
end

