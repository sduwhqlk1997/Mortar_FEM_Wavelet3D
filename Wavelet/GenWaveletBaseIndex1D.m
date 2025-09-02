function [BaseFunIndex,IntegralBlock] = GenWaveletBaseIndex1D(interval,j0,j,BaseType)
%GENWAVELETBASEINDEX1D 生成一维小波基的信息表和各基函数对应的积分分段（保证在每个子区间上均为多项式）
%   BaseFunIndex：一维小波基函数信息表
%       1：SpaceType，0表示尺度基，1表示小波基
%       2：level
%       3：基函数编号
%       4：支集左端点
%       5：支集右端点
%       6：基函数坐标索引
%       7：基函数编号
%   IntegralBlock：一维小波基函数的积分分段，为矩阵，每一行表示一个基函数，非NaN数据表示节点（从左到右依次排列）
%       顺序与BaseFunIndex一致

switch BaseType
    case "linear"

    case "quadratic"
        % 各生成函数的支集
        supp_phi = [0,3];
        supp_phi_b1 = [0,1];
        supp_phi_b2 = [0,2];
        supp_psi = [0,3];
        supp_psi_bD = [0,2.5];
        % 各生成函数的分段节点
        %   尺度基
        NodePhi=[0,1,2,3]; % 内部基
        NodePhiB1=[0,1]; % 左边界基1
        NodePhiB2=[0,1,2]; % 左边界基2
        %   小波基
        NodePsi=[0,0.5,1,1.5,2,2.5,3]; % 内部基
        NodePsiB=[0,0.5,1,1.5,2,2.5]; % 左边界基

        BaseFunIndex=zeros(2^j+2,7);
        IntegralBlock=NaN(2^j+2,7);

        BaseFunIndex(1:2^j0+2,1:6) = [0,j0,1,supp_phi_b1.*2^(-j0),0;...
            0,j0,2,supp_phi_b2.*2^(-j0),2/3*2^(-j0);...
            zeros(2^j0-2,1),j0.*ones(2^j0-2,1),(3:2^j0)',(supp_phi+(0:2^j0-3)').*2^(-j0),2^(-j0)*((3:2^j0)'-2);...
            0,j0,2^j0+1,1-supp_phi_b2(end:-1:1)*(2^-j0),1-2/3*2^(-j0);...
            0,j0,2^j0+2,1-supp_phi_b1(end:-1:1)*(2^-j0),1];
%[1-2^(-j0+1),1]    [1-2^(-j0),1]
        IntegralBlock(1,1:2)=NodePhiB1.*2^(-j0);
        IntegralBlock(2,1:3)=NodePhiB2.*2^(-j0);
        IntegralBlock(3:2^j0,1:4)=(NodePhi+(0:2^j0-3)').*2^(-j0);
        IntegralBlock(2^j0+1,1:3)=1-NodePhiB2(end:-1:1)*2^(-j0);
        IntegralBlock(2^j0+2,1:2)=1-NodePhiB1(end:-1:1)*2^(-j0);

        if j>=j0+1
            for i = j0:j-1
                BaseFunIndex(2^(i)+3:2^(i+1)+2,1:6) = [1,i,1,supp_psi_bD.*2^(-i),2.^(-i-1);...
                    ones(2^i-2,1),i.*ones(2^i-2,1),(2:2^i-1)',(supp_psi+(0:2^i-3)').*2^(-i),2.^(-i-1)+2.^(-i).*((2:2^i-1)'-1);...
                    1,i,2^i,[1-5/2^(i+1),1],1-2.^(-i-1)];

                IntegralBlock(2^(i)+3,1:6)=NodePsiB.*2^(-i);
                IntegralBlock(2^(i)+4:2^(i+1)+1,1:7)=(NodePsi+(0:2^i-3)').*2^(-i);
                IntegralBlock(2^(i+1)+2,1:6)=1-NodePsiB(end:-1:1).*2^(-i);
            end
        end
    otherwise
        error("基函数类型未定义")
end
BaseFunIndex(:,4:6) = (interval(2)-interval(1)).*BaseFunIndex(:,4:6)+interval(1);
BaseFunIndex(:,end) = 1:size(BaseFunIndex,1); % 给每个基函数的编号
IntegralBlock=(interval(2)-interval(1)).*IntegralBlock+interval(1);
end

