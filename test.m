%% 测试 NNZIndexFEMLag = FindNNPosFEMLag(FEMRegionMeshIndex2D,LagMeshIndex)
clear
clc
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet\FEM_FEM';
addpath(genpath(currentPath));
baseType="linear";
oriFEM=[0;0;0];
a=2;b=2;c=2;
NxFEM=5; NyFEM=5; NzFEM=3;
intervalLag=[0.5,1.5;0.5,1.5];
zFixed=2;
NxLag=5; NyLag=5;
x_step=linspace(intervalLag(1,1),intervalLag(1,2),NxLag).';
y_step=linspace(intervalLag(2,1),intervalLag(2,2),NyLag).';

[P,T,Pb,Tb] = genMesh3D(oriFEM,a,b,c,NxFEM,NyFEM,NzFEM,baseType);
viewMesh(P,T)
Myfindnode(Pb)
[RegionBaseIndex2D,RegionMeshIndex2D] = Extract2DRegion(Pb,Tb,[intervalLag;zFixed,zFixed]);
[LagBaseIndex,LagMeshIndex] = GenLagBaseIndex(baseType,x_step,y_step);
NNZIndexFEMLag = FindNNPosFEMLag(RegionMeshIndex2D,LagMeshIndex);