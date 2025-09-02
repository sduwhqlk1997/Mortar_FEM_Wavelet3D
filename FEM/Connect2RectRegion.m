function [K,F,Dof_index] = Connect2RectRegion(K1,K2,Dof_index1,Dof_index2,F1,F2)
%CONNECT2RECTREGION 拼接两个矩形区域
%   Dof_Index：第一列自由度类型1:u,2:v,3:w,4:phi；2-4列坐标
ax=norm([Dof_index1(:,2);Dof_index2(:,2)],'inf');
ay=norm([Dof_index1(:,3);Dof_index2(:,3)],'inf');
az=norm([Dof_index1(:,4);Dof_index2(:,4)],'inf');
[interface1,interface2]=find(Dof_index1(:,1)==Dof_index2(:,1).'&...
    abs(Dof_index1(:,2)-Dof_index2(:,2).')/ax<1e-10&...
    abs(Dof_index1(:,3)-Dof_index2(:,3).')/ay<1e-10&...
    abs(Dof_index1(:,4)-Dof_index2(:,4).')/az<1e-10);
interface2=interface2+size(K1,1);
K=blkdiag(K1,K2);
Dof_index=[Dof_index1;Dof_index2];

K(interface1,:)=K(interface1,:)+K(interface2,:);
K(:,interface1)=K(:,interface1)+K(:,interface2);

K(interface2,:)=[];
K(:,interface2)=[];
Dof_index(interface2,:)=[];

if exist('F1','var')&&exist('F2','var')
    F=[F1;F2];
    F(interface1,:)=F(interface1,:)+F(interface2,:);
    F(interface2,:)=[];
else
    F=[];
end
end

