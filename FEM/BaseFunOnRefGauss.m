function [phi,Dphi] = BaseFunOnRefGauss(order,FEMBaseType)
%������������䵼���ڲο���Ԫ[-1,1]^3�ϵ�Gauss�㴦��ֵ
%   orderΪһ��ά����Gauss��ĸ���
%   FEMBaseType������������
%   phi �洢���л�������Gauss�㴦��ֵ��ÿһ�ж�Ӧһ��������
%   Dphi Ϊһ��num*order^3*3����ά���飬������ά�ȵ�1��2��3��ֱ��Ӧ�Ź���x,y,z�ĵ�����ÿһ��Ķ�Ӧ��ϵ��phi��ͬ
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

