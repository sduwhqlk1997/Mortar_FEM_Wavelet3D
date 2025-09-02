function y = GaussIntForWavelet(f,interval_Sol,interval,order,dim,BaseIndex_x,DofX_num,BaseIndex_y,DofY_num,BaseIndex_z,DofZ_num)
%GAUSSINTFORWAVELET ר�ż��㱻��������С�����Ļ���
%   ����˵��ͬGenGaussBlock.m��GaussIntegralBlocks.m
if dim==1
    X = GenWaveGaussBlock(interval_Sol,interval,dim,BaseIndex_x,DofX_num);
    y = GaussIntegralBlocks(f,order,dim,X);
elseif dim==2
    [X,Y] = GenWaveGaussBlock(interval_Sol,interval,dim,BaseIndex_x,DofX_num,BaseIndex_y,DofY_num);
    y = GaussIntegralBlocks(f,order,dim,X,Y);
elseif dim==3
    [X,Y,Z] = GenWaveGaussBlock(interval_Sol,interval,dim,BaseIndex_x,DofX_num,BaseIndex_y,DofY_num,BaseIndex_z,DofZ_num);
    y = GaussIntegralBlocks(f,order,dim,X,Y,Z);
else
    error("ά��ҪС��3")
end
end

