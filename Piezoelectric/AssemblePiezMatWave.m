function [K,M,BaseIndex,DofIndex] = AssemblePiezMatWave(Mat1DRef,BaseFunIndex_ref1D,...
    jcoef,J,j0,MaterialPara,interval)
%ASSEMBLEPIEZMATWAVE 用小波基组装压电基底区域的刚度矩阵和载荷向量
% 输入：
%   Mat1DRef：[0,1]区间上的一维参考矩阵，为cell数组。1：K00_1D，2：K10_1D，3：K11_1D
%   BaseFunIndex：[0,1]区间上的一维小波基函数信息表
%       1：SpaceType，0表示尺度基，1表示小波基
%       2：level
%       3：基函数编号
%       4：支集左端点
%       5：支集右端点
%       6：基函数坐标索引
%       7：基函数编号
%   Jcoef,J：j1coef*j1+j2coef*j2+j3coef*j3<=J
%   j0：最低的level
%   MaterialPara：cell数组，材料参数
%       1. 弹性张量 2. 压电张量 3. 介电常数 4. 密度
%   interval：器件区域，为3x2矩阵，每行表示一个维度的区间
% 输出：
%   K：刚度矩阵，M：质量矩阵
%   BaseIndex：基函数信息表
%       1：type2_x      2：type2_y   3: type2_z （尺度基or小波基）;
%       4：j_x      5：j_y   6: j_z
%       7：k_x      8：k_y   9: k_z
%       10~11：interval_x     12~13：interval_y  14~15: interval_z
%       16~18：坐标索引
%       19~24: 该基函数的支集，依次为x,y,z方向的基函数
%   DofIndex：自由度索引，自由度与基函数的对应关系。
%       1：自由度与基函数对应关系（与BaseIndex的行对应）
%       2：基函数类型
% 生成三维稀疏基函数信息
[DofIndex_sparse,Table_sparse] = GenWaveletSparseInfo3D(j0,J,...
    BaseFunIndex_ref1D,BaseFunIndex_ref1D,BaseFunIndex_ref1D,jcoef);
% 生成三维参考基本矩阵
K00_1D=cell2mat(Mat1DRef(1)); K10_1D=cell2mat(Mat1DRef(2)); K11_1D=cell2mat(Mat1DRef(3));
A00=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K00_1D,K00_1D,K00_1D);
A11=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K11_1D,K00_1D,K00_1D);
A12=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K10_1D,K10_1D.',K00_1D);
A22=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K00_1D,K11_1D,K00_1D);
A13=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K10_1D,K00_1D,K10_1D.');
A23=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K00_1D,K10_1D,K10_1D.');
A33=AssembleWaveMatrix3D(DofIndex_sparse,Table_sparse,K00_1D,K00_1D,K11_1D);
% 通过仿射变换形成总体基本矩阵
JacobiDet=prod(interval(:,2)-interval(:,1));
JacobiInv=1./(interval(:,2)-interval(:,1));
% 将基本总体矩阵变换到指定区域
A00=JacobiDet*A00;
A11=JacobiDet*JacobiInv(1)^2*A11;
A12=JacobiDet*JacobiInv(1)*JacobiInv(2)*A12;
A22=JacobiDet*JacobiInv(2)^2*A22;
A13=JacobiDet*JacobiInv(1)*JacobiInv(3)*A13;
A23=JacobiDet*JacobiInv(2)*JacobiInv(3)*A23;
A33=JacobiDet*JacobiInv(3)^2*A33;
% 
A12=A12.';A13=A13.';A23=A23.';
% 提取材料参数
c_p=cell2mat(MaterialPara(1));
e_p=cell2mat(MaterialPara(2));
epcl_p=cell2mat(MaterialPara(3));
rho_p=cell2mat(MaterialPara(4));
% 组装刚度矩阵和质量矩阵
N=size(A00,1);
M = rho_p*A00;
M = blkdiag(M,M,M,sparse(N,N));
% Kuu
u11 = C(c_p,1,1,1,1)*A11+C(c_p,1,1,1,2)*(A12+A12.')+C(c_p,1,1,1,3)*(A13+A13.')+ ...
    C(c_p,1,2,1,2)*A22+C(c_p,1,2,1,3)*(A23+A23.')+ ...
    C(c_p,1,3,1,3)*A33;
u12 = C(c_p,1,1,1,2)*A11+C(c_p,1,1,2,2)*A12.'+C(c_p,1,1,2,3)*A13.'+ ...
    C(c_p,1,2,1,2)*A12+C(c_p,1,2,2,2)*A22+C(c_p,1,2,2,3)*A23.'+ ...
    C(c_p,1,3,1,2)*A13+C(c_p,1,3,2,2)*A23+C(c_p,1,3,2,3)*A33;
u13 = C(c_p,1,1,1,3)*A11+C(c_p,1,1,2,3)*A12.'+C(c_p,1,1,3,3)*A13.'+ ...
    C(c_p,1,2,1,3)*A12+C(c_p,1,2,2,3)*A22+C(c_p,1,2,3,3)*A23.'+ ...
    C(c_p,1,3,1,3)*A13+C(c_p,1,3,2,3)*A23+C(c_p,1,3,3,3)*A33;
u22 = C(c_p,1,2,1,2)*A11+C(c_p,1,2,2,2)*(A12+A12.')+C(c_p,1,2,2,3)*(A13+A13.')+ ...
    C(c_p,2,2,2,2)*A22+C(c_p,2,2,2,3)*(A23+A23.')+ ...
    C(c_p,2,3,2,3)*A33;
u23 = C(c_p,1,2,1,3)*A11+C(c_p,1,2,2,3)*A12.'+C(c_p,1,2,3,3)*A13.'+ ...
    C(c_p,2,2,1,3)*A12+C(c_p,2,2,2,3)*A22+C(c_p,2,2,3,3)*A23.'+ ...
    C(c_p,2,3,1,3)*A13+C(c_p,2,3,2,3)*A23+C(c_p,2,3,3,3)*A33;
u33 = C(c_p,1,3,1,3)*A11+C(c_p,1,3,2,3)*(A12+A12.')+C(c_p,1,3,3,3)*(A13+A13.')+ ...
    C(c_p,2,3,2,3)*A22+C(c_p,2,3,3,3)*(A23+A23.')+ ...
    C(c_p,3,3,3,3)*A33;
Kuu = [u11,u12,u13;u12.',u22,u23;u13.',u23.',u33];
% Kup
phi11 = E(e_p,1,1,1)*A11+E(e_p,2,1,1)*A12.'+E(e_p,3,1,1)*A13.'+...
    E(e_p,1,2,1)*A12+E(e_p,2,2,1)*A22+E(e_p,3,2,1)*A23.'+...
    E(e_p,1,3,1)*A13+E(e_p,2,3,1)*A23+E(e_p,3,3,1)*A33;
phi21 = E(e_p,1,1,2)*A11+E(e_p,2,1,2)*A12.'+E(e_p,3,1,2)*A13.'+...
    E(e_p,1,2,2)*A12+E(e_p,2,2,2)*A22+E(e_p,3,2,2)*A23.'+...
    E(e_p,1,3,2)*A13+E(e_p,2,3,2)*A23+E(e_p,3,3,2)*A33;
phi31 = E(e_p,1,1,3)*A11+E(e_p,2,1,3)*A12.'+E(e_p,3,1,3)*A13.'+...
    E(e_p,1,2,3)*A12+E(e_p,2,2,3)*A22+E(e_p,3,2,3)*A23.'+...
    E(e_p,1,3,3)*A13+E(e_p,2,3,3)*A23+E(e_p,3,3,3)*A33;
Kup = [phi11;phi21;phi31];
% Kp
Kp = epcl_p(1,1)*A11+epcl_p(1,2)*A12+epcl_p(1,3)*A13+...
    epcl_p(2,1)*A12.'+epcl_p(2,2)*A22+epcl_p(2,3)*A23+...
    epcl_p(3,1)*A13.'+epcl_p(3,2)*A23.'+epcl_p(3,3)*A33;
% 总体刚度矩阵
K=[Kuu,Kup;...
    Kup.',-Kp];
% 生成自由度索引表
BaseFunIndex_x=BaseFunIndex_ref1D;
BaseFunIndex_y=BaseFunIndex_ref1D;
BaseFunIndex_z=BaseFunIndex_ref1D;

BaseFunIndex_x(:,4:6)=BaseFunIndex_x(:,4:6)*(interval(1,2)-interval(1,1))+interval(1,1);
BaseFunIndex_y(:,4:6)=BaseFunIndex_y(:,4:6)*(interval(2,2)-interval(2,1))+interval(2,1);
BaseFunIndex_z(:,4:6)=BaseFunIndex_z(:,4:6)*(interval(3,2)-interval(3,1))+interval(3,1);

BaseIndex=[BaseFunIndex_x(DofIndex_sparse(:,1),1),... type2
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
    BaseFunIndex_x(DofIndex_sparse(:,1),4:5),... support set
    BaseFunIndex_y(DofIndex_sparse(:,2),4:5),...
    BaseFunIndex_z(DofIndex_sparse(:,3),4:5)];
DofIndex=(1:size(BaseIndex,1))';
DofIndex=[repmat(DofIndex,4,1),kron([1;2;3;4],ones(N,1))];
end

