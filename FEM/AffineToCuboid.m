function [Jacobi,Jacobi_inv,y] = AffineToCuboid(x,P,T,flag)
%���ο���Ԫ[-1,1]^3�ĵ�ӳ�䵽�����嵥ԪT��
%   xΪ�ο���Ԫ[-1,1]^3�ϵ�ĳ����
%   PΪ�����ȫ������
%   TΪ��Ԫ�ڵ��������Ķ�Ӧ��ϵ
%   ���yΪһ��3ά���飬������ά�ȵ�ÿһ��ֱ��Ӧ�任����x,y,z����,ÿһ���ÿһ�ж�Ӧһ����Ԫ
%   Jacobi:����Ԫ����任��Ӧ��Jacobi����ʽ,ÿ�ж�Ӧһ����Ԫ
%   Jacobi_inv:ÿһ�ж�Ӧһ��Ԫ�أ���һ��Ϊdref_x/dx,�ڶ���Ϊdref_y/dy,������Ϊdref_z/dz
%   ���flagΪ0����Ҫ��xӳ�䵽���е�Ԫ��
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

