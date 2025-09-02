function IntegralBlock = GenWaveletIntInterval1D(interval,BaseType,SpaceType,j,k)
%GENWAVELETINTINTERVAL1D 按小波的分段生成积分区间的节点，保证在每个子区间上均为多项式（一维小波）
%   interval：定义域
%   BaseType：基函数类型
%   SpaceType：空间类型
%   j：Level
%   k：基函数编号

switch BaseType
    case "linear"
    case "quadratic"
        switch SpaceType
            case 0 % 尺度基
                NodePhi=[0;1;2;3]; % 内部基
                NodePhiB1=[0;1]; % 左边界基1
                NodePhiB1=[0;1;2]; % 左边界基2
                if k==1

                elseif k==2

                elseif k>=3 && k<=2^j

                elseif k==2^j+1

                elseif k==2^j+2
                    
                else
                    error("k超出范围,该基函数不存在")
                end
            case 1 % 小波基
                NodePsi=[0;0.5;1;1.5;2;2.5;3]; % 内部基
                NodePsiB=[0;0.5;1;1.5;2;2.5]; % 左边界基
        end
end
end

