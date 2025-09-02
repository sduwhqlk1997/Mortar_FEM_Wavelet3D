clear
%% 测试 Extract2DRegion
clc
clear
currentPath = 'D:\Code\M\Mortar_FEM_Wavelet\FEM_FEM\FEM';
addpath(genpath(currentPath));
ori=[0;0;0];
a=2; b=2; c=2;
Nx=5; Ny=5; Nz=3;
type="linear";
[P,T,Pb,Tb] = genMesh3D(ori,a,b,c,Nx,Ny,Nz,type);
viewMesh(P,T)
Myfindnode(Pb)
[RegionBaseIndex2D,RegionMeshIndex2D] = Extract2DRegion(Pb,Tb,[0.2,1.8;0.2,1.8;2,2]);
%% 生成线性基函数的索引矩阵
clc;
clear;
ori=[1;1;1];
a=1;b=1;c=1;
type="linear";
Nx=2;Ny=2;Nz=2;
[P,T,Pb,Tb] = genMesh3D(ori,a,b,c,Nx,Ny,Nz,type);
viewMesh(P,T)
Myfindnode(Pb(Tb,:))
% figure(2)
% viewMesh(P,T)
% Myfindnode(Pb)
%% 测试y = FEMBaseFunRef(type,x,index,dx,dy,dz)
clear;
clc;
[X,Y,Z]=meshgrid([-1;1]);
X=X(:); Y=Y(:); Z=Z(:);
for i = 1:8
    value=FEMBaseFunRef("linear",X,Y,Z,i,0,0,0);
end
%% 测试 value = FEMBaseFun(X,baseType,ele,baseNum,diff)
clear;
clc;
[X,Y,Z]=meshgrid([0;2],[1;4],[2;6]);
X=X(:); Y=Y(:); Z=Z(:);
for i = 1:8
    % value=FEMBaseFunRef("linear",[X,Y,Z],i,0,0,0);
    value=FEMBaseFun(X,Y,Z,"linear",[0,2;1,4;2,6],i,[0;0;0]);
    1
end
%% 测试基函数的和
x=-1:0.1:1;
[x,y,z]=meshgrid(x);
f=@(x,y,z,index)FEMBaseFunRef("quadratic",x,y,z,index,0,0,0);
ff=zeros(size(x,1),size(x,2));
for i=1:27
    ff=ff+f(x,y,z,i);
end