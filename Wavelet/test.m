%% 测试y = WaveletBaseRef(BaseType,SpaceType,x,j,k,diff)
clear
clc
BaseType="quadratic";
SpaceType=0;
x=0:0.001:1;
j=2;
diff=1;
figure(1)
hold on
title("尺度基")
for k=1:2^j+2
    y = WaveletBaseRef(BaseType,SpaceType,x,j,k,diff);
    plot(x,y)
end
hold off
SpaceType=1;
j=4;
diff=1;
figure(2)
hold on
title("小波基")
for k=1:2^j
    y = WaveletBaseRef(BaseType,SpaceType,x,j,k,diff);
    plot(x,y)
end
hold off
%% 测试BaseFunIndex = GenWaveletBaseIndex1D(interval,j0,j,BaseType)
clear
BaseType="quadratic";
j0=2;
j=3;
interval=[3,4];
[BaseFunIndex,IntegralBlock] = GenWaveletBaseIndex1D(interval,j0,j,BaseType);
%% 测试WaveletRegionBaseIndex2D = Extract2DRegionWavelet(WaveIndex4Lag,interval)
clear
BaseType="quadratic";
J=6;
j0=2;
j=J-2*j0;
interval=[1,2;3,4;5,6];
[BaseFunIndex_x,IntegralBlock_x] = GenWaveletBaseIndex1D(interval(1,:),j0,j,BaseType);
[BaseFunIndex_y,IntegralBlock_y] = GenWaveletBaseIndex1D(interval(2,:),j0,j,BaseType);
[BaseFunIndex_z,IntegralBlock_z] = GenWaveletBaseIndex1D(interval(3,:),j0,j,BaseType);
[DofIndex_sparse,Table_sparse] = GenWaveletSparseInfo3D(j0,J,...
    BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z);
WaveIndex4Lag = GenWaveBaseIndex4Lag(DofIndex_sparse,...
    BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z);
block=[1.25,1.75;3.25,3.75;6,6];
WaveletRegionBaseIndex2D = Extract2DRegionWavelet(WaveIndex4Lag,block);
suppSelect=[BaseFunIndex_x(DofIndex_sparse(WaveletRegionBaseIndex2D,1),4:5),...
    BaseFunIndex_y(DofIndex_sparse(WaveletRegionBaseIndex2D,2),4:5),...
    BaseFunIndex_z(DofIndex_sparse(WaveletRegionBaseIndex2D,3),4:5)];
index=(1:size(DofIndex_sparse,1))';
index(WaveletRegionBaseIndex2D)=[];
suppDel=[BaseFunIndex_x(DofIndex_sparse(index,1),4:5),...
    BaseFunIndex_y(DofIndex_sparse(index,2),4:5),...
    BaseFunIndex_z(DofIndex_sparse(index,3),4:5)];