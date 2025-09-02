function K = AssembleFEMMatrix(coeff,Pb,Tb,trail,test,phi,Dphi,Jacobi,Jacobi_inv,w)
%ASSEMBLEFEMMATRIX 此处显示有关此函数的摘要
%   此处显示详细说明
%组装{int(coeff*f*g)}刚度矩阵
%   
%   Pb,Tb分别为有限元节点和网格信息矩阵
%   trail和test均为三维向量，分别储存试探函数和测试函数关于x,y,z的导数阶
%   order为高斯积分一个维度上的Gauss点数
%   coeff为系数，若为变系数要输入所有Gauss点的系数值，每一行代表一个单元，每一列代表一个Gauss点，若为常系数只要输入系数的标量值即可
%   phi和Dphi分别为参考基函数和其导数在Gauss节点的值，这里假设被积函数至多为一阶导数,且trail和test要么全为0阶导数要么全为某个变量的一阶导数
%   Jacobi,Jacobi_inv分别为各单元的Jacobi行列式和Jacobi矩阵的逆矩阵
%   w为Gauss系数
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
%             if i~=j                                                        %考虑对称性
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

