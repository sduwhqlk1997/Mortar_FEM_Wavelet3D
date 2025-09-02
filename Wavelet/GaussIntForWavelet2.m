function y = GaussIntForWavelet2(f,order,dim,IntegralBlock1_x,IntegralBlock2_x,...
    IntegralBlock1_y,IntegralBlock2_y,...
    IntegralBlock1_z,IntegralBlock2_z)
%GAUSSINTFORWAVELET2 专门计算被积函数含小波基的积分
%   此处显示详细说明
switch dim
    case 1
        X = GenWaveGaussBlock2(dim,IntegralBlock1_x,IntegralBlock2_x);
        y = GaussIntegralBlocks(f,order,dim,X);
    case 2
        [X,Y] = GenWaveGaussBlock2(dim,IntegralBlock1_x,IntegralBlock2_x,...
            IntegralBlock1_y,IntegralBlock2_y);
        y = GaussIntegralBlocks(f,order,dim,X,Y);
    case 3
        [X,Y,Z] = GenWaveGaussBlock2(dim,IntegralBlock1_x,IntegralBlock2_x,...
            IntegralBlock1_y,IntegralBlock2_y,...
            IntegralBlock1_z,IntegralBlock2_z);
        y = GaussIntegralBlocks(f,order,dim,X,Y,Z);
end
end

