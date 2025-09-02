function [P,T,Pb,Tb] = genMesh2D(ori,a,b,Nx,Ny,type,x,y)
%GENMESH2D 二维矩形区域的均匀网格剖分
%   ori:原点位置,a:沿x轴正向的长度,b:沿y轴正向的长度
%   Nx: 沿x方向的剖分节点数,Ny:沿y方向的剖分节点数
%   type: 剖分类型，9:张量积形式的双二次元网格
if nargin == 6
    x = linspace(ori(1),ori(1)+a,Nx);
    y = linspace(ori(2),ori(2)+b,Ny);
end
[X,Y] = meshgrid(x,y);
P = [X(:),Y(:)];
index = 1:Nx*Ny;
index = reshape(index,Ny,Nx);
T = [reshape(index(1:Ny-1,1:Nx-1),[],1),reshape(index(2:Ny,1:Nx-1),[],1),reshape(index(2:Ny,2:Nx),[],1),reshape(index(1:Ny-1,2:Nx),[],1)];
if type==9
    P1 = (P(reshape(index(1:Ny,1:Nx-1),[],1),:)+P(reshape(index(1:Ny,2:Nx),[],1),:))/2;
    P2 = (P(reshape(index(1:Ny-1,1:Nx),[],1),:)+P(reshape(index(2:Ny,1:Nx),[],1),:))/2;
    P3 = (P(reshape(index(1:Ny-1,1:Nx-1),[],1),:)+P(reshape(index(2:Ny,2:Nx),[],1),:))/2;
    index2 = (size(P,1)+1):(size(P,1)+size(P1,1));
    index2 = reshape(index2,Ny,[]);
    Tb = [T,reshape(index2(1:Ny-1,1:Nx-1),[],1),reshape(index2(2:Ny,1:Nx-1),[],1)];
    Pb = [P;P1];
    index3 = (size(Pb,1)+1):(size(Pb,1)+size(P2,1));
    index3 = reshape(index3,[],Nx);
    Tb = [Tb,reshape(index3(1:Ny-1,1:Nx-1),[],1),reshape(index3(1:Ny-1,2:Nx),[],1)];
    Pb = [Pb;P2];
    index4 = (size(Pb,1)+1):(size(Pb,1)+size(P3,1));
%     index4 = reshape(index4,[],Nx-1);
    Tb = [Tb,index4'];
    Pb = [Pb;P3];
end
end

