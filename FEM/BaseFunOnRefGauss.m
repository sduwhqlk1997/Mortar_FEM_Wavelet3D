function [phi,Dphi] = BaseFunOnRefGauss(order,FEMBaseType)
%计算基函数及其导数在参考单元[-1,1]^3上的Gauss点处的值
%   order为一个维度上Gauss点的个数
%   FEMBaseType：基函数类型
%   phi 存储所有基函数在Gauss点处的值，每一行对应一个基函数
%   Dphi 为一个num*order^3*3的三维数组，第三个维度的1，2，3块分别对应着关于x,y,z的导数，每一块的对应关系和phi相同
switch FEMBaseType
    case "linear"
        num=8;
    case "quadratic"
        num=27;
end
[pt,w] = genRefGauss3DCube(order);
phi = zeros(num,length(w));
DphiDx = phi;
DphiDy = phi;
DphiDz = phi;
for i=1:num
    phi(i,:) = FEMBaseFunRef(FEMBaseType,pt(:,1),pt(:,2),pt(:,3),i,0,0,0)';
    DphiDx(i,:) = FEMBaseFunRef(FEMBaseType,pt(:,1),pt(:,2),pt(:,3),i,1,0,0)';
    DphiDy(i,:) = FEMBaseFunRef(FEMBaseType,pt(:,1),pt(:,2),pt(:,3),i,0,1,0)';
    DphiDz(i,:) = FEMBaseFunRef(FEMBaseType,pt(:,1),pt(:,2),pt(:,3),i,0,0,1)';
end
Dphi(:,:,1) =  DphiDx;
Dphi(:,:,2) =  DphiDy;
Dphi(:,:,3) =  DphiDz;
end

