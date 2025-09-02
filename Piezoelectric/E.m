function eijk = E(coef,i,j,k)
%E 在压缩存储后的压电张量中提取参数
%   coef：压缩存储后的压电张量
%   ijk:压缩存储前参数对应的下标
if j==k
    m=j;
elseif j+k==5
    m=4;
elseif j+k==4
    m=5;
elseif j+k==3
    m=6;
end
eijk = coef(i,m);
end

