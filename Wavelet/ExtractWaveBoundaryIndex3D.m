function boundary_index = ExtractWaveBoundaryIndex3D(domain_B,Dof_Index)
%ExtractBoundaryIndex3D 根据选定的边界区域，提取相应的自由度索引
%   boundary_index：边界区域包含的自由度索引，与Dof_Index相对应
%   domain_B：边界区域，为3*2矩阵，第一行表示x方向，第二行表示y方向，第三行表示z方向
%             若某一方向左右端点值相等，表示该边界在该方向的坐标恒为该值
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

index=(1:size(Dof_Index,1))';
if domain_B(1,1)==domain_B(1,2) % 该边界平行于yz面
    boundary_index=index(Dof_Index(:,16)==domain_B(1,1) &...
        Dof_Index(:,17)>=domain_B(2,1) & Dof_Index(:,17)<=domain_B(2,2)&...
        Dof_Index(:,18)>=domain_B(3,1) & Dof_Index(:,18)<=domain_B(3,2));
elseif domain_B(2,1)==domain_B(2,2) % 该边界平行于xz面
    boundary_index=index(Dof_Index(:,17)==domain_B(2,1) &...
        Dof_Index(:,16)>=domain_B(1,1) & Dof_Index(:,16)<=domain_B(1,2)&...
        Dof_Index(:,18)>=domain_B(3,1) & Dof_Index(:,18)<=domain_B(3,2));
elseif domain_B(3,1)==domain_B(3,2) % 该边界平行于xy面
    boundary_index=index(Dof_Index(:,18)==domain_B(3,1) &...
        Dof_Index(:,16)>=domain_B(1,1) & Dof_Index(:,16)<=domain_B(1,2)&...
        Dof_Index(:,17)>=domain_B(2,1) & Dof_Index(:,17)<=domain_B(2,2));
end
end

