function [X,Y,Z] = GenWaveGaussBlock2(dim,IntegralBlock1_x,IntegralBlock2_x,...
    IntegralBlock1_y,IntegralBlock2_y,...
    IntegralBlock1_z,IntegralBlock2_z)
%GENWAVEGAUSSBLOCK2 生成积分区间段
%   此处显示详细说明
if dim ==1 || dim ==2 || dim ==3
    IntegralBlock1_x=IntegralBlock1_x(:);
    IntegralBlock2_x=IntegralBlock2_x(:);
    X=intersect2Sequence(IntegralBlock1_x,IntegralBlock2_x);
end
if dim ==2 || dim ==3
    IntegralBlock1_y=IntegralBlock1_y(:);
    IntegralBlock2_y=IntegralBlock2_y(:);
    Y=intersect2Sequence(IntegralBlock1_y,IntegralBlock2_y);
end
if dim ==3
    IntegralBlock1_z=IntegralBlock1_z(:);
    IntegralBlock2_z=IntegralBlock2_z(:);
    Z=intersect2Sequence(IntegralBlock1_z,IntegralBlock2_z);
end
    function sequence= intersect2Sequence(sequence1,sequence2) % 合并两个积分序列的相交部分，此时两个序列中都不含NaN
        interval=[max([sequence1(1),sequence2(1)]),min([sequence1(end),sequence2(end)])];
        Block=unique([sequence1;sequence2]);
        sequence=ExtractBlock(interval,Block);
        function IntBlock=ExtractBlock(interval,Block) % 截取Block中落在interval区间中的区间段
            left=Block<=interval(1);
            right=Block>=interval(2);
            ileft=any(left); % 若区间左端点不在任何区间段内则为0
            iright=any(right); % 若区间右端点不在任何区间段内则为0
            Block(left|right)=[]; % 删掉不在区间中的元素
            if  ileft&&iright  % 区间在区间段内
                IntBlock=[interval(1);Block;interval(2)];
            elseif ~ileft&&iright
                IntBlock=[Block;interval(2)];
            elseif ileft&&~iright
                IntBlock=[interval(1);Block];
            elseif ~ileft&&~iright
                IntBlock=Block;
            end
        end
    end
end

