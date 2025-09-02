function F = AssembleFEMVector(f,Pb,Tb,phi,Jacobi,y,w)
%组装载荷向量(\int(f*phi))
%   f(x,y,z)为右端项
%   Pb,Tb为有限元节点和有限元网格信息
%   phi为Gauss点处的基函数值,这里假设基函数各变量的导数阶均为0
%   Jacobi为各单元Jacobi行列式
%   w为参考单元上的Gauss权重
%   y各单元上Gauss点的坐标
%[Jacobi,Jacobi_inv,y] = AffineToCuboid(x,P,T,flag)
F1 = reshape(f(reshape(y(:,:,1),[],1),reshape(y(:,:,2),[],1),reshape(y(:,:,3),[],1)),size(Tb,1),[]);%计算f在所有单元上的Gauss点的值,每行对应一个单元
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

