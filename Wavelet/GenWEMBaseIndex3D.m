function WaveletDofIndex3D = GenWEMBaseIndex3D(interval,DofIndex_sparse,...
    BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z)
%GENWAVELETBASEINDEX3D 生成三维小波基的信息表
%   WaveletDofIndex3D：从左到右依次为:
        %1：type2_x      2：type2_y   3: type2_z （尺度基or小波基）;
        %4：j_x      5：j_y   6: j_z
        %7：k_x      8：k_y   9: k_z
        %10~11：interval_x     12~13：interval_y  14~15: interval_z
        %16~18：坐标索引
        %19~20：与该基函数x方向做过Matching的1D基函数的定义区间
        %21~22：与该基函数y方向做过Matching的1D基函数的定义区间 
        %23~24：与该基函数z方向做过Matching的1D基函数的定义区间 (若两个位置全为0则该方向未做过Matching)
        %25~30: 该基函数的支集，依次为x,y,z方向的基函数
%   interval：该组基函数的定义区域，第一行表示x方向，第二行表示y方向，第三行表示z方向
%   DofIndex_sparse:稀疏基索引，1~3列分别为对应x,y,z方向一维小波基的编号，与BaseFunIndex对应
%   BaseFunIndex_x,BaseFunIndex_y,BaseFunIndex_z:分别为x,y,z方向的一维参考基的基函数信息，
%       1：SpaceType，0表示尺度基，1表示小波基
%       2：level
%       3：基函数编号
%       4：支集左端点
%       5：支集右端点
%       6：基函数坐标索引
%       7：基函数编号

WaveletDofIndex3D=[BaseFunIndex_x(DofIndex_sparse(:,1),1),... type2
    BaseFunIndex_y(DofIndex_sparse(:,2),1),...
    BaseFunIndex_z(DofIndex_sparse(:,3),1),... 
    BaseFunIndex_x(DofIndex_sparse(:,1),2),... j
    BaseFunIndex_y(DofIndex_sparse(:,2),2),...
    BaseFunIndex_z(DofIndex_sparse(:,3),2),...
    BaseFunIndex_x(DofIndex_sparse(:,1),3),... k
    BaseFunIndex_y(DofIndex_sparse(:,2),3),...
    BaseFunIndex_z(DofIndex_sparse(:,3),3),...
    kron(ones(size(DofIndex_sparse,1),1),interval(1,:)),... interval
    kron(ones(size(DofIndex_sparse,1),1),interval(2,:)),...
    kron(ones(size(DofIndex_sparse,1),1),interval(3,:)),...
    BaseFunIndex_x(DofIndex_sparse(:,1),end-1),... coordinate index
    BaseFunIndex_y(DofIndex_sparse(:,2),end-1),...
    BaseFunIndex_z(DofIndex_sparse(:,3),end-1),...
    zeros(size(DofIndex_sparse,1),6),... Matching information
    BaseFunIndex_x(DofIndex_sparse(:,1),4:5),... support set
    BaseFunIndex_y(DofIndex_sparse(:,2),4:5),...
    BaseFunIndex_z(DofIndex_sparse(:,3),4:5)];

end

