function K = AssembleFEMMatrix(coeff,Pb,Tb,trail,test,phi,Dphi,Jacobi,Jacobi_inv,w)
%ASSEMBLEFEMMATRIX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%��װ{int(coeff*f*g)}�նȾ���
%   
%   Pb,Tb�ֱ�Ϊ����Ԫ�ڵ��������Ϣ����
%   trail��test��Ϊ��ά�������ֱ𴢴���̽�����Ͳ��Ժ�������x,y,z�ĵ�����
%   orderΪ��˹����һ��ά���ϵ�Gauss����
%   coeffΪϵ������Ϊ��ϵ��Ҫ��������Gauss���ϵ��ֵ��ÿһ�д���һ����Ԫ��ÿһ�д���һ��Gauss�㣬��Ϊ��ϵ��ֻҪ����ϵ���ı���ֵ����
%   phi��Dphi�ֱ�Ϊ�ο����������䵼����Gauss�ڵ��ֵ��������豻����������Ϊһ�׵���,��trail��testҪôȫΪ0�׵���ҪôȫΪĳ��������һ�׵���
%   Jacobi,Jacobi_inv�ֱ�Ϊ����Ԫ��Jacobi����ʽ��Jacobi����������
%   wΪGaussϵ��
ii = zeros(size(Tb,1),size(Tb,2)*size(Tb,2));
jj = zeros(size(Tb,1),size(Tb,2)*size(Tb,2));
kk = zeros(size(Tb,1),size(Tb,2)*size(Tb,2));
state = 1;
if sum(trail)==0&&sum(test)==0
    for i = 1:size(Tb,2)
        for j = 1:size(Tb,2)
            kk(:,state) = coeff.*(phi(i,:).*phi(j,:))*w.*Jacobi;
%             ii(:,state) = Tb(:,i);
%             jj(:,state) = Tb(:,j);
            jj(:,state) = Tb(:,i);
            ii(:,state) = Tb(:,j);
            state = state+1;
%             if i~=j                                                        %���ǶԳ���
%                 kk(:,state) = kk(:,state-1);
%                 ii(:,state) = jj(:,state-1);
%                 jj(:,state) = ii(:,state-1);
%                 state = state+1;
%             end
        end
    end
else
    def_trail=find(trail==1);
    def_test=find(test==1);
    for i = 1:size(Tb,2)
        for j = 1:size(Tb,2)
            kk(:,state) = coeff.*Dphi(i,:,def_trail).*Dphi(j,:,def_test)*w.*Jacobi.*Jacobi_inv(:,def_trail).*Jacobi_inv(:,def_test);
%             ii(:,state) = Tb(:,i);
%             jj(:,state) = Tb(:,j);
            jj(:,state) = Tb(:,i);
            ii(:,state) = Tb(:,j);
            state = state+1;
        end
    end
end
K = sparse(ii(:),jj(:),kk(:),size(Pb,1),size(Pb,1));
% K(abs(K)<1e-20)=0;
end

