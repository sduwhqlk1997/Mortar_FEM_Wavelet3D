function f = ApproxWaveFun3D(x,y,z,u,Dof_Index,diff,type1,DofType)
%ApproxWaveFun3D 数值解
%   x,y,z：要计算的点的三个方向坐标
%   u：基函数的系数，顺序与Dof_Index对应
%   Dof_Index：从左到右依次为:
%       1：type2_x      2：type2_y   3: type2_z （尺度基or小波基）;
%       4：j_x      5：j_y   6: j_z
%       7：k_x      8：k_y   9: k_z
%       10~11：interval_x     12~13：interval_y  14~15: interval_z
%       16~18：坐标索引
%       19~24: 该基函数的支集，依次为x,y,z方向的基函数
%       25：自由度类型
%   diff: 导数阶
%   DofType：自由度类型
if exist('DofType','var')
    Dof_Index(Dof_Index(:,25)~=DofType,:)=[];
end
N = size(Dof_Index,1);
f=zeros(size(x,1),size(x,2));
for i = 1:N
    f=f+...
        u(i)*WaveletBaseFun3D(x,y,z,...
        [Dof_Index(i,10:11);Dof_Index(i,12:13);Dof_Index(i,14:15)],...
        Dof_Index(i,4:6),Dof_Index(i,7:9),diff,type1,Dof_Index(i,1:3));
end
end

