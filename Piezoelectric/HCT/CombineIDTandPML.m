function [K,Dof_indexHCT] ...
    = CombineIDTandPML(K,Dof_indexHCT,K_L,K_R,Dof_indexHCT_L,Dof_indexHCT_R)
%COMBINEIDTANDPML 将级联后的IDT矩阵与PML的schur补矩阵拼接
%   K,Dof_indexHCT：IDT的索引和矩阵
%   K_L,K_R,Dof_indexHCT_L,Dof_indexHCT_R：左（L）、右（R）PML区域的索引和矩阵

% 左
% ay=norm([Dof_indexHCT(:,2);Dof_indexHCT_L(:,2)],'inf');
% az=norm([Dof_indexHCT(:,3);Dof_indexHCT_L(:,3)],'inf');
ay=1e-7;
az=1e-6;
[i,j]=find(abs(Dof_indexHCT(:,2)-Dof_indexHCT_L(:,2).')/ay<1e-10&...
    abs(Dof_indexHCT(:,3)-Dof_indexHCT_L(:,3).')/az<1e-10&...
    Dof_indexHCT(:,4)==Dof_indexHCT_L(:,4).'&...
    Dof_indexHCT(:,1)==-1);
K(i,i)=K(i,i)+K_L(j,j);
% 右
[i,j]=find(abs(Dof_indexHCT(:,2)-Dof_indexHCT_R(:,2).')/ay<1e-10&...
    abs(Dof_indexHCT(:,3)-Dof_indexHCT_R(:,3).')/az<1e-10&...
    Dof_indexHCT(:,4)==Dof_indexHCT_R(:,4).'&...
    Dof_indexHCT(:,1)==-2);
K(i,i)=K(i,i)+K_R(j,j);
end

