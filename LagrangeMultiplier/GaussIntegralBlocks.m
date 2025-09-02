function y = GaussIntegralBlocks(f,order,dim,X,Y,Z)
%GAUSSINTEGRALBLOCKS 分段Gauss积分
%   f:被积函数
%   order:Gauss积分的阶
%   dim:积分的维度
%   X,Y,Z三个方向的分点,均为列向量
if order==4
    x = [0.861136311594053;-0.861136311594053;0.339981043584856;-0.339981043584856];
    h = [0.347854845137454;0.347854845137454;0.652145154862546;0.652145154862546];
    
elseif order==8
    h=[0.1012285363,0.1012285363,0.2223810345,0.2223810345,0.3137066459,0.3137066459,0.3626837834,0.3626837834]';
    x=[0.9602898565,-0.9602898565,0.7966664774,-0.7966664774,0.5255324099,-0.5255324099,0.1834346425,-0.1834346425]';
elseif order==2
    h = [1;1];
    x = [1/sqrt(3);-1/sqrt(3)];
end
if dim==1
    interval = [X(1:end-1),X(2:end)];
    a = (interval(:,2)-interval(:,1))/2;
    mid = (interval(:,2)+interval(:,1))/2;
    x=kron(x,a)+kron(ones(size(x,1),1),mid);
    x=reshape(x,size(interval,1),[]);
    x=x.';
    y = sum(f(x).'*h.*a);
elseif dim==2
    % 生成Gauss点
    [X_G,Y_G]=meshgrid(x,x);
    x = [X_G(:),Y_G(:)];
    [h1,h2] = meshgrid(h,h);
    w = reshape(h1.*h2,[],1);
    % 生成积分分片
    Nx=length(X);Ny=length(Y);
    [X_mesh,Y_mesh] = meshgrid(X,Y);
    P = [X_mesh(:),Y_mesh(:)];
    index = 1:Nx*Ny;
    index = reshape(index,Ny,Nx);
    T = [reshape(index(1:Ny-1,1:Nx-1),[],1),reshape(index(2:Ny,1:Nx-1),[],1),reshape(index(2:Ny,2:Nx),[],1),reshape(index(1:Ny-1,2:Nx),[],1)];
    mid = (P(T(:,1),:)+P(T(:,3),:))./2;
    a = abs(P(T(:,1),1)-P(T(:,4),1))./2;
    b = abs(P(T(:,1),2)-P(T(:,2),2))./2;
    xx = reshape(kron(x(:,1),a)+kron(ones(size(x,1),1),mid(:,1)),size(T,1),[]);
    yy = reshape(kron(x(:,2),b)+kron(ones(size(x,1),1),mid(:,2)),size(T,1),[]);
    y = sum(f(xx,yy)*w.*a.*b);
elseif dim==3
    % 生成Gauss点
    [X_G,Y_G,Z_G]=meshgrid(x);
    x = [X_G(:),Y_G(:),Z_G(:)];
    [h1,h2,h3] = meshgrid(h);
    w = reshape(h1.*h2.*h3,[],1);
    % 生成积分分片
    Nx=length(X);Ny=length(Y);Nz = length(Z);
    [X_mesh,Y_mesh,Z_mesh] = meshgrid(X,Y,Z);
    P = [X_mesh(:),Y_mesh(:),Z_mesh(:)];
    index = 1:Nx*Ny*Nz;
    index = reshape(index,Ny,Nx,Nz);
    T = [reshape(index(1:Ny-1,1:Nx-1,1:Nz-1),[],1),...
        reshape(index(2:Ny,1:Nx-1,1:Nz-1),[],1),...
        reshape(index(1:Ny-1,2:Nx,1:Nz-1),[],1),...
        reshape(index(2:Ny,2:Nx,1:Nz-1),[],1),...
        reshape(index(1:Ny-1,1:Nx-1,2:Nz),[],1),...
        reshape(index(2:Ny,1:Nx-1,2:Nz),[],1),...
        reshape(index(1:Ny-1,2:Nx,2:Nz),[],1),...
        reshape(index(2:Ny,2:Nx,2:Nz),[],1)];
    mid = (P(T(:,1),:)+P(T(:,8),:))./2;
    a = abs(P(T(:,1),1)-P(T(:,3),1))./2;
    b = abs(P(T(:,1),2)-P(T(:,2),2))./2;
    c = abs(P(T(:,1),3)-P(T(:,5),3))./2;
    xx = reshape(kron(x(:,1),a)+kron(ones(size(x,1),1),mid(:,1)),size(T,1),[]);
    yy = reshape(kron(x(:,2),b)+kron(ones(size(x,1),1),mid(:,2)),size(T,1),[]);
    zz = reshape(kron(x(:,3),c)+kron(ones(size(x,1),1),mid(:,3)),size(T,1),[]);
    y = sum(f(xx,yy,zz)*w.*a.*b.*c);
else
    error("维数至多为三维")
end
end

