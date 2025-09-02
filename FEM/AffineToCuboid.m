function [Jacobi,Jacobi_inv,y] = AffineToCuboid(x,P,T,flag)
%将参考单元[-1,1]^3的点映射到长方体单元T上
%   x为参考单元[-1,1]^3上的某个点
%   P为网格点全局坐标
%   T为单元节点与网格点的对应关系
%   输出y为一个3维数组，第三个维度的每一层分别对应变换后点的x,y,z坐标,每一层的每一行对应一个单元
%   Jacobi:各单元仿射变换对应的Jacobi行列式,每行对应一个单元
%   Jacobi_inv:每一行对应一个元素，第一列为dref_x/dx,第二列为dref_y/dy,第三列为dref_z/dz
%   如果flag为0则不需要把x映射到所有单元上
mid = (P(T(:,1),:)+P(T(:,8),:))./2;
a = abs(P(T(:,1),1)-P(T(:,3),1))./2;
b = abs(P(T(:,1),2)-P(T(:,2),2))./2;
c = abs(P(T(:,1),3)-P(T(:,5),3))./2;
if flag == 1
    y(:,:,1) = reshape(kron(x(:,1),a)+kron(ones(size(x,1),1),mid(:,1)),size(T,1),[]);
    y(:,:,2) = reshape(kron(x(:,2),b)+kron(ones(size(x,1),1),mid(:,2)),size(T,1),[]);
    y(:,:,3) = reshape(kron(x(:,3),c)+kron(ones(size(x,1),1),mid(:,3)),size(T,1),[]);
end
Jacobi = a.*b.*c;
Jacobi_inv = [1./a,1./b,1./c];
end

