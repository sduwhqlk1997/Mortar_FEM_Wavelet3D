function [P,T,Pb,Tb] = genMesh3D(ori,a,b,c,Nx,Ny,Nz,type,x,y,z)
%生成长方体区域的二次网格
%   ori为长宽高分别平行于x,y,z轴的长方体的原点，a,b,c分别为长宽高,Nx,Ny,Nz,分别为沿x,y,轴划分的格点数
%   输出P，T为网格格点全局坐标和网格单元局部坐标对应的全局编号
%   Pb,Tb分别为二次元的格点坐标和有限元网格格点编号
%   type为要使用基函数的类型，20为20点二次元，27为27点二次元
%   x,y,z为三个方向的网格步长向量，若省略则默认进行均匀剖分
if nargin == 8
    x = linspace(ori(1),ori(1)+a,Nx);
    y = linspace(ori(2),ori(2)+b,Ny);
    z = linspace(ori(3),ori(3)+c,Nz);
end
[X,Y,Z] = meshgrid(x,y,z);
if isempty(Nx)&& isempty(Ny) &&isempty(Nz)
    Nx=length(x); Ny=length(y); Nz=length(z);
end
P = [X(:),Y(:),Z(:)];
index = 1:Nx*Ny*Nz;
index = reshape(index,Ny,Nx,Nz);
T = [reshape(index(1:Ny-1,1:Nx-1,1:Nz-1),[],1),reshape(index(2:Ny,1:Nx-1,1:Nz-1),[],1),reshape(index(1:Ny-1,2:Nx,1:Nz-1),[],1),reshape(index(2:Ny,2:Nx,1:Nz-1),[],1),reshape(index(1:Ny-1,1:Nx-1,2:Nz),[],1),reshape(index(2:Ny,1:Nx-1,2:Nz),[],1),reshape(index(1:Ny-1,2:Nx,2:Nz),[],1),reshape(index(2:Ny,2:Nx,2:Nz),[],1)];
switch type
    case "linear"
        Pb=P;
        Tb=T;
    case "quadratic"
        % 生成各边中点
        P1 = (P(reshape(index(:,1:Nx-1,:),[],1),:)+P(reshape(index(:,2:Nx,:),[],1),:))/2;
        P2 = (P(reshape(index(1:Ny-1,:,:),[],1),:)+P(reshape(index(2:Ny,:,:),[],1),:))/2;
        P3 = (P(reshape(index(:,:,1:Nz-1),[],1),:)+P(reshape(index(:,:,2:Nz),[],1),:))/2;
        % 生成面中点
        P4 = (P(reshape(index(1:Ny,1:Nx-1,1:Nz-1),[],1),:)+P(reshape(index(1:Ny,2:Nx,2:Nz),[],1),:))/2;
        P5 = (P(reshape(index(1:Ny-1,1:Nx,1:Nz-1),[],1),:)+P(reshape(index(2:Ny,1:Nx,2:Nz),[],1),:))/2;
        P6 = (P(reshape(index(1:Ny-1,1:Nx-1,1:Nz),[],1),:)+P(reshape(index(2:Ny,2:Nx,1:Nz),[],1),:))/2;
        % 生成体中点
        P7 = (P(reshape(index(1:Ny-1,1:Nx-1,1:Nz-1),[],1),:)+P(reshape(index(2:Ny,2:Nx,2:Nz),[],1),:))/2;
        Pb = [P;P1;P2;P3;P4;P5;P6;P7];
        len = size(P,1);
        index1 = len+1:len+size(P1,1); index1 = reshape(index1,Ny,Nx-1,Nz);%x方向边
        len  = len+size(P1,1);
        index2 = len+1:len+size(P2,1); index2 = reshape(index2,Ny-1,Nx,Nz);%y方向边
        len = len+size(P2,1);
        index3 = len+1:len+size(P3,1); index3 = reshape(index3,Ny,Nx,Nz-1);%z方向边
        len = len+size(P3,1);
        index4 = len+1:len+size(P4,1); index4 = reshape(index4,Ny,Nx-1,Nz-1);   %左右面中点
        len=len+size(P4,1);
        index5 = len+1:len+size(P5,1); index5 = reshape(index5,Ny-1,Nx,Nz-1);   %后前面中点
        len = len+size(P5,1);
        index6 = len+1:len+size(P6,1); index6 = reshape(index6,Ny-1,Nx-1,Nz);   %下上面中点
        len = len+size(P6,1);
        index7 = len+1:len+size(P7,1);%体中点
        Tb = [T,reshape(index1(1:Ny-1,1:Nx-1,1:Nz-1),[],1),reshape(index1(2:Ny,1:Nx-1,1:Nz-1),[],1),reshape(index1(1:Ny-1,1:Nx-1,2:Nz),[],1),reshape(index1(2:Ny,1:Nx-1,2:Nz),[],1),...
            reshape(index2(1:Ny-1,1:Nx-1,1:Nz-1),[],1),reshape(index2(1:Ny-1,2:Nx,1:Nz-1),[],1),reshape(index2(1:Ny-1,1:Nx-1,2:Nz),[],1),reshape(index2(1:Ny-1,2:Nx,2:Nz),[],1),...
            reshape(index3(1:Ny-1,1:Nx-1,1:Nz-1),[],1),reshape(index3(1:Ny-1,2:Nx,1:Nz-1),[],1),reshape(index3(2:Ny,1:Nx-1,1:Nz-1),[],1),reshape(index3(2:Ny,2:Nx,1:Nz-1),[],1), ...
            reshape(index4(1:Ny-1,1:Nx-1,1:Nz-1),[],1),reshape(index4(2:Ny,1:Nx-1,1:Nz-1),[],1),...
            reshape(index5(1:Ny-1,1:Nx-1,1:Nz-1),[],1),reshape(index5(1:Ny-1,2:Nx,1:Nz-1),[],1),...
            reshape(index6(1:Ny-1,1:Nx-1,1:Nz-1),[],1),reshape(index6(1:Ny-1,1:Nx-1,2:Nz),[],1),...
            index7'];
end
end

