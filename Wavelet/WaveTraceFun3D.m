function [u_B,index_I,index_B] = WaveTraceFun3D(f,boundary_index,BaseIndex,type1)
%TRACEFUN3D 处理某一未知量的Dirichlet边界条件
%   f：Dirichlet边界条件
%   K：刚度矩阵
%   F：载荷向量
%   u_B：边界自由度的近似解
%   index_I: 非边界自由度索引，编号与Dof_Index的行号对应
%   index_B: 边界自由度索引，编号与Dof_Index的行号对应
%   boundary_index：自由度编号，与Dof_Index的行坐标对应
%   BaseIndex：从左到右依次为:
%       1：type2_x      2：type2_y   3: type2_z （尺度基or小波基）;
%       4：j_x      5：j_y   6: j_z
%       7：k_x      8：k_y   9: k_z
%       10~11：interval_x     12~13：interval_y  14~15: interval_z
%       16~18：坐标索引
%       19~24: 该基函数的支集，依次为x,y,z方向的基函数

index_I=(1:size(BaseIndex,1))';
A = zeros(size(boundary_index,1));

for i = 1:size(boundary_index,1)
    A(:,i)=WaveletBaseFun3D(BaseIndex(boundary_index,16),...
        BaseIndex(boundary_index,17),BaseIndex(boundary_index,18),...
        [BaseIndex(boundary_index(i),10:11);BaseIndex(boundary_index(i),12:13);BaseIndex(boundary_index(i),14:15)],...
        BaseIndex(boundary_index(i),4:6),BaseIndex(boundary_index(i),7:9),...
        [0,0,0],type1,BaseIndex(boundary_index(i),1:3));
    % WaveBaseFunMatching3D(Dof_Index(boundary_index,16),...
    %     Dof_Index(boundary_index,17),Dof_Index(boundary_index,18),...
    %     [Dof_Index(boundary_index(i),10:11);Dof_Index(boundary_index(i),12:13);Dof_Index(boundary_index(i),14:15)],...
    %     Dof_Index(boundary_index(i),4:6),Dof_Index(boundary_index(i),7:9),...
    %     [0,0,0],type1,Dof_Index(boundary_index(i),1:3),...
    %     [Dof_Index(boundary_index(i),19:20);...
    %     Dof_Index(boundary_index(i),21:22);...
    %     Dof_Index(boundary_index(i),23:24)]);
end
F = f(BaseIndex(boundary_index,16),BaseIndex(boundary_index,17),BaseIndex(boundary_index,18));
u_B = A\F;
index_I(boundary_index,:)=[];
index_B = boundary_index;
end

