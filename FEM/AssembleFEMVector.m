function F = AssembleFEMVector(f,Pb,Tb,phi,Jacobi,y,w)
%��װ�غ�����(\int(f*phi))
%   f(x,y,z)Ϊ�Ҷ���
%   Pb,TbΪ����Ԫ�ڵ������Ԫ������Ϣ
%   phiΪGauss�㴦�Ļ�����ֵ,�������������������ĵ����׾�Ϊ0
%   JacobiΪ����ԪJacobi����ʽ
%   wΪ�ο���Ԫ�ϵ�GaussȨ��
%   y����Ԫ��Gauss�������
%[Jacobi,Jacobi_inv,y] = AffineToCuboid(x,P,T,flag)
F1 = reshape(f(reshape(y(:,:,1),[],1),reshape(y(:,:,2),[],1),reshape(y(:,:,3),[],1)),size(Tb,1),[]);%����f�����е�Ԫ�ϵ�Gauss���ֵ,ÿ�ж�Ӧһ����Ԫ
ii = zeros(size(Tb,1),size(Tb,2));
ff_1=ii;
state = 1;
for i = 1:size(Tb,2)
    ff_1(:,state) = F1.*phi(i,:)*w.*Jacobi;
    ii(:,state) = Tb(:,i);
    state = state+1;
end
F = accumarray(ii(:),ff_1(:),[size(Pb,1),1]);
end

