function [K,Dof_index] = MatchingTwoCubeWaveDomain(K1,K2,Dof_index1,Dof_index2,num_unknows)
%MatchingTwoCubeDomain 拼接两个不相交的长方体区域 (只能用于单个面的matching)
%   Dof_Index：从左到右依次为:
        %1：type2_x      2：type2_y   3: type2_z;
        %4：j_x      5：j_y   6: j_z
        %7：k_x      8：k_y   9: k_z
        %10~11：interval_x     12~13：interval_y  14~15: interval_z
        %16~18：坐标索引
        %19~20：与该基函数x方向做过Matching的1D基函数的定义区间
        %21~22：与该基函数y方向做过Matching的1D基函数的定义区间 
        %23~24：与该基函数z方向做过Matching的1D基函数的定义区间 (若两个位置全为0则该方向未做过Matching)
        %25~30: 该基函数的支集，以此为x,y,z方向的基函数
%   num_unknows: 方程未知量个数
%   Dof_index1,Dof_index2: 分别为K1和K2的自由度信息矩阵
if ~exist("num_unknows",'var')
    num_unknows=1;
end

N1 = size(Dof_index1,1); N2 = size(Dof_index2,1); N=N1+N2;
[interface1,interface2]=find((abs(Dof_index1(:,16)-Dof_index2(:,16).')./norm(Dof_index1(:,16),inf))<1e-10 &...
    (abs(Dof_index1(:,17)-Dof_index2(:,17).')./norm(Dof_index1(:,17),inf))<1e-10 &...
    (abs(Dof_index1(:,18)-Dof_index2(:,18).')./norm(Dof_index1(:,18),inf))<1e-10);
interface_p=Dof_index1(interface1,16:18);

if interface_p(1,1)==interface_p(end,1) % 该交界面垂直于x轴
    flag=19;
elseif interface_p(1,2)==interface_p(end,2) % 该交界面平行于y轴
    flag=21;
elseif interface_p(1,3)==interface_p(end,3) % 该交界面平行于z轴
    flag=23;
end
K = blkdiag(K1,K2);
if num_unknows>=2 % 当变量数大于1个时要调整矩阵排列，使得矩阵按不同自由度依次排列
    reorder=zeros(1,size(K,1));
    for i = 1:num_unknows % 生成矩阵重排向量
        reorder(N*(i-1)+1:N*i)=[(i-1)*N1+1:i*N1,size(K1,1)+((i-1)*N2+1):(size(K1,1)+i*N2)];
    end
    K = K(reorder,reorder);
end
Dof_index=[Dof_index1;Dof_index2];
Dof_index(interface1,flag:flag+1)=Dof_index2(interface2,flag-9:flag-8); % 标记做过的matching操作
Dof_index(interface1,flag+6:flag+7)=[min([Dof_index(interface1,flag+6),Dof_index2(interface2,flag+6)],[],2),...
    max([Dof_index(interface1,flag+7),Dof_index2(interface2,flag+7)],[],2)]; % Matching后基函数的支集也要做相应变换

reserve_index=reshape(interface1+(0:num_unknows-1)*N,[],1); % 需要保留的行列
delete_index=reshape(N1+interface2+(0:num_unknows-1)*N,[],1); % 需要删除的行列

K(reserve_index,:) = K(reserve_index,:) + K(delete_index,:);
K(:,reserve_index) = K(:,reserve_index) + K(:,delete_index);
K(delete_index,:)=[];
K(:,delete_index)=[];
Dof_index(N1+interface2,:) = [];
end

