function value = WaveBaseFunMatching3D(x,y,z,...
    interval,j,k,diff,type1,type2,Matching)
%BASEFUNMATCHING3D Matching后的基函数
%   x,y,z：要计算函数值的点的x,y,z坐标
%   interval：计算区域，第一行表示x区间，第二行表示y区间，第三行表示z区间
%   j：level,为3*1矩阵，分别为j_x，j_y，j_z
%   k：3*1矩阵，分别为k_x,k_y,k_z
%   diff：导数阶，3*1矩阵，分别表示x、y和z的导数阶，最高为各1阶
%   type1：基函数类型，为字符串
%   type2：尺度基(0)or小波基(1)，为3*1矩阵
%   Matching：3*2向量，与该基函数Matching过的基函数的定义区间，
%       第一行表示x方向，第二行表示y方向，第二行表示z方向

if type1=="linear"
    k_end = @(type2,j) (type2==0).*(2^j+1)+(type2==1).*2^j;
elseif type1=="quadratic"
    k_end = @(type2,j) (type2==0).*(2^j+2)+(type2==1).*2^j;
end
f_x =@(x) WaveletBaseFun1D(x,interval(1,:),j(1),k(1),diff(1),type1,type2(1));
f_y =@(y) WaveletBaseFun1D(y,interval(2,:),j(2),k(2),diff(2),type1,type2(2));
f_z =@(z) WaveletBaseFun1D(z,interval(3,:),j(3),k(3),diff(3),type1,type2(3));
if Matching(1,1)<Matching(1,2) % x方向做了Matching
    
    if interval(1,1)==Matching(1,2) % 如果是左端点和右端点Matching
        f_x=@(x) f_x(x).*(x>=interval(1,1)&x<=interval(1,2))+...
            WaveletBaseFun1D(x,Matching(1,:),j(1),k_end(type2(1),j(1)),diff(1),type1,type2(1)).*(x>=Matching(1,1)&x<Matching(1,2));
    elseif interval(1,2)==Matching(1,1) % 如果是右端点和左端点Matching
        f_x=@(x) f_x(x).*(x>=interval(1,1)&x<=interval(1,2))+...
            WaveletBaseFun1D(x,Matching(1,:),j(1),1,diff(1),type1,type2(1)).*(x>Matching(1,1)&x<=Matching(1,2));
    end
end

if Matching(2,1)<Matching(2,2) % y方向做了Matching
    
    if interval(2,1)==Matching(2,2) % 如果是左端点和右端点Matching
        f_y=@(y) f_y(y).*(y>=interval(2,1)&y<=interval(2,2))+...
            WaveletBaseFun1D(y,Matching(2,:),j(2),k_end(type2(2),j(2)),diff(2),type1,type2(2)).*(y>=Matching(2,1)&y<Matching(2,2));
    elseif interval(2,2)==Matching(2,1) % 如果是右端点和左端点Matching
        f_y=@(y) f_y(y).*(y>=interval(2,1)&y<=interval(2,2))+...
            WaveletBaseFun1D(y,Matching(2,:),j(2),1,diff(2),type1,type2(2)).*(y>Matching(2,1)&y<=Matching(2,2));
    end
end

if Matching(3,1)<Matching(3,2) % y方向做了Matching
    
    if interval(3,1)==Matching(3,2) % 如果是左端点和右端点Matching
        f_z=@(z) f_z(z).*(z>=interval(3,1)&z<=interval(3,2))+...
            WaveletBaseFun1D(z,Matching(3,:),j(3),k_end(type2(3),j(3)),diff(3),type1,type2(3)).*(z>=Matching(3,1)&z<Matching(3,2));
    elseif interval(3,2)==Matching(3,1) % 如果是右端点和左端点Matching
        f_z=@(z) f_z(z).*(z>=interval(3,1)&z<=interval(3,2))+...
            WaveletBaseFun1D(z,Matching(3,:),j(3),1,diff(3),type1,type2(3)).*(z>Matching(3,1)&z<=Matching(3,2));
    end
end
value = f_x(x).*f_y(y).*f_z(z);
end

