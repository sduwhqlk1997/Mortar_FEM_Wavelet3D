function F = AssembleWaveVectorMatching3D(f,Dof_index,type_test,order)
%AssembleVectorMatching3D 用Matching后的基函数组装载荷向量
%   f：右端项
%   Dof_Index：从左到右依次为:
%1：type2_x      2：type2_y   3: type2_z;
%4：j_x      5：j_y   6: j_z
%7：k_x      8：k_y   9: k_z
%10~11：interval_x     12~13：interval_y  14~15: interval_z
%16~18：坐标索引
%19~20：与该基函数x方向做过Matching的1D基函数的定义区间
%21~22：与该基函数y方向做过Matching的1D基函数的定义区间
%23~24：与该基函数z方向做过Matching的1D基函数的定义区间 (若两个位置全为0则该方向未做过Matching)
%25~30: 该基函数的支集，以此为x,y,z方向的基函数
%   type_test：基函数类型
%   order：Gauss函数阶数

N = size(Dof_index,1);
F = zeros(N,1);

for i = 1:N
    % i
    fun=@(x,y,z) f(x,y,z).*...
        WaveletBaseFun1D(x,Dof_index(i,10:11),Dof_index(i,4),Dof_index(i,7),0,type_test,Dof_index(i,1)).*...
        WaveletBaseFun1D(y,Dof_index(i,12:13),Dof_index(i,5),Dof_index(i,8),0,type_test,Dof_index(i,2)).*...
        WaveletBaseFun1D(z,Dof_index(i,14:15),Dof_index(i,6),Dof_index(i,9),0,type_test,Dof_index(i,3));
    supp_x=GenSuppofWavelet1D(Dof_index(i,10:11),type_test,Dof_index(i,1),Dof_index(i,4),Dof_index(i,7));
    supp_y=GenSuppofWavelet1D(Dof_index(i,12:13),type_test,Dof_index(i,2),Dof_index(i,5),Dof_index(i,8));
    supp_z=GenSuppofWavelet1D(Dof_index(i,14:15),type_test,Dof_index(i,3),Dof_index(i,6),Dof_index(i,9));
    F(i) = GaussIntForWavelet(fun,...
        [Dof_index(i,10:11);Dof_index(i,12:13);Dof_index(i,14:15)],...
        [supp_x;supp_y;supp_z],...
        order,3,[0,Dof_index(i,4)],1,[0,Dof_index(i,5)],1,[0,Dof_index(i,6)],1);
    if type_test=="linear"
        k_end = @(type2,j,interval1,interval2) (interval1(1)==interval2(2))...
            .*((type2==0).*(2^j+1)+(type2==1).*2^j)+(interval1(2)==interval2(1)).*1;
    elseif type_test=="quadratic"
        k_end = @(type2,j,interval1,interval2) (interval1(1)==interval2(2))...
            .*((type2==0).*(2^j+2)+(type2==1).*2^j)+(interval1(2)==interval2(1)).*1;
    end
    if Dof_index(i,19)<Dof_index(i,20) % 只有x方向做了matching
  
        k_x=k_end(Dof_index(i,1),Dof_index(i,4),Dof_index(i,10:11),Dof_index(i,19:20));
        % +--
        fun=@(x,y,z) f(x,y,z).*...
            WaveletBaseFun1D(x,Dof_index(i,19:20),Dof_index(i,4),k_x,0,type_test,Dof_index(i,1)).*...
            WaveletBaseFun1D(y,Dof_index(i,12:13),Dof_index(i,5),Dof_index(i,8),0,type_test,Dof_index(i,2)).*...
            WaveletBaseFun1D(z,Dof_index(i,14:15),Dof_index(i,6),Dof_index(i,9),0,type_test,Dof_index(i,3));

        supp_x=GenSuppofWavelet1D(Dof_index(i,19:20),type_test,Dof_index(i,1),Dof_index(i,4),k_x);
        supp_y=GenSuppofWavelet1D(Dof_index(i,12:13),type_test,Dof_index(i,2),Dof_index(i,5),Dof_index(i,8));
        supp_z=GenSuppofWavelet1D(Dof_index(i,14:15),type_test,Dof_index(i,3),Dof_index(i,6),Dof_index(i,9));

        F(i) = F(i)+GaussIntForWavelet(fun,...
            [Dof_index(i,19:20);Dof_index(i,12:13);Dof_index(i,14:15)],...
            [supp_x;supp_y;supp_z],...
            order,3,[0,Dof_index(i,4)],1,[0,Dof_index(i,5)],1,[0,Dof_index(i,6)],1);
    elseif Dof_index(i,21)<Dof_index(i,22) % 只有y方向做了matching
        
        k_y=k_end(Dof_index(i,2),Dof_index(i,5),Dof_index(i,12:13),Dof_index(i,21:22));
        % -+-
        fun=@(x,y,z) f(x,y,z).*...
            WaveletBaseFun1D(x,Dof_index(i,10:11),Dof_index(i,4),Dof_index(i,7),0,type_test,Dof_index(i,1)).*...
            WaveletBaseFun1D(y,Dof_index(i,21:22),Dof_index(i,5),k_y,0,type_test,Dof_index(i,2)).*...
            WaveletBaseFun1D(z,Dof_index(i,14:15),Dof_index(i,6),Dof_index(i,9),0,type_test,Dof_index(i,3));

        supp_x=GenSuppofWavelet1D(Dof_index(i,10:11),type_test,Dof_index(i,1),Dof_index(i,4),Dof_index(i,7));
        supp_y=GenSuppofWavelet1D(Dof_index(i,21:22),type_test,Dof_index(i,2),Dof_index(i,5),k_y);
        supp_z=GenSuppofWavelet1D(Dof_index(i,14:15),type_test,Dof_index(i,3),Dof_index(i,6),Dof_index(i,9));

        F(i) = F(i)+GaussIntForWavelet(fun,...
            [Dof_index(i,10:11);Dof_index(i,21:22);Dof_index(i,14:15)],...
            [supp_x;supp_y;supp_z],...
            order,3,[0,Dof_index(i,4)],1,[0,Dof_index(i,5)],1,[0,Dof_index(i,6)],1);
    elseif Dof_index(i,23)<Dof_index(i,24) % 只有z方向做了matching
        
        k_z=k_end(Dof_index(i,3),Dof_index(i,6),Dof_index(i,14:15),Dof_index(i,23:24));
        % --+
        fun=@(x,y,z) f(x,y,z).*...
            WaveletBaseFun1D(x,Dof_index(i,10:11),Dof_index(i,4),Dof_index(i,7),0,type_test,Dof_index(i,1)).*...
            WaveletBaseFun1D(y,Dof_index(i,12:13),Dof_index(i,5),Dof_index(i,8),0,type_test,Dof_index(i,2)).*...
            WaveletBaseFun1D(z,Dof_index(i,23:24),Dof_index(i,6),k_z,0,type_test,Dof_index(i,3));

        supp_x=GenSuppofWavelet1D(Dof_index(i,10:11),type_test,Dof_index(i,1),Dof_index(i,4),Dof_index(i,7));
        supp_y=GenSuppofWavelet1D(Dof_index(i,12:13),type_test,Dof_index(i,2),Dof_index(i,5),Dof_index(i,8));
        supp_z=GenSuppofWavelet1D(Dof_index(i,23:24),type_test,Dof_index(i,3),Dof_index(i,6),k_z);

        F(i) = F(i)+GaussIntForWavelet(fun,...
            [Dof_index(i,10:11);Dof_index(i,12:13);Dof_index(i,23:24)],...
            [supp_x;supp_y;supp_z],...
            order,3,[0,Dof_index(i,4)],1,[0,Dof_index(i,5)],1,[0,Dof_index(i,6)],1);
    end
end
end

