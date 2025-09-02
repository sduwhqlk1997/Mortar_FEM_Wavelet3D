function [K,F,u_index_I,v_index_I,w_index_I,phi_index_I,u_B,v_B,w_B,phi_B] ...
    = TreatPiezDiriBoundWave3D(K,F,Diri_info,BaseIndex,WaveBaseType)
%TREATPIEZDIRIBOUNDWAVE3D 处理压电方程的Dirichlet边界条件
%   K,F：刚度矩阵和载荷向量
%   u_index_I,v_index_I,w_index_I,phi_index_I：n*2矩阵，
%       1：为该自由度对应的刚度矩阵中的位置
%       2：为对应的Dof_Index中的位置
%   u_B,v_B,w_B,phi_B：边界上自由度的近似解
%       1：该自由度对应的Dof_Index中的位置
%       2：为该自由度的近似解
%   Diri_info：cell数组，每一行表示一组边界条件信息，
%       1：为对应的自由度编号(1:u 2:v 3:w 4:phi)
%       2：为边界点索引集合，与Dof_Index相对应
%       3：为边界上的精确解函数 (齐次边界直接设置为0即可)
%       4: 0: 齐次边界 1：非齐次边界
%   BaseIndex：基函数信息表
%       1：type2_x      2：type2_y   3: type2_z （尺度基or小波基）;
%       4：j_x      5：j_y   6: j_z
%       7：k_x      8：k_y   9: k_z
%       10~11：interval_x     12~13：interval_y  14~15: interval_z
%       16~18：坐标索引
%       19~24: 该基函数的支集，依次为x,y,z方向的基函数
%   WaveBaseType：基函数类型

u_B = [];   v_B = [];   w_B = [];   phi_B = [];

for i=1:size(Diri_info,1)
    switch cell2mat(Diri_info(i,4))
        case 1 % 非齐次边界
            [u_b,~,index_b] = ...
                WaveTraceFun3D(cell2mat(Diri_info(i,3)),cell2mat(Diri_info(i,2)),BaseIndex,WaveBaseType);
        case 0 % 齐次边界
            u_b = zeros(size(cell2mat(Diri_info(i,2)),1),1);
            index_b = cell2mat(Diri_info(i,2));

    end
    switch cell2mat(Diri_info(i,1))
        case 1 % u
            u_B = [u_B;[index_b,u_b]];
        case 2 % v
            v_B = [v_B;[index_b,u_b]];
        case 3 % v
            w_B = [w_B;[index_b,u_b]];
        case 4 % phi
            phi_B = [phi_B;[index_b,u_b]];
        otherwise
            error('该未知量不存在')
    end
end
N = size(BaseIndex,1);
index_delete = [];
u_index_I = (1:N)'; v_index_I = u_index_I; w_index_I = u_index_I; phi_index_I = u_index_I;
if ~isempty(u_B)
    F = F-sum(u_B(:,2).'.*K(:,u_B(:,1)),2);
    index_delete = [index_delete;u_B(:,1)];
    u_index_I(u_B(:,1),:) = [];
    u_index_I = [(1:size(u_index_I,1))',u_index_I];
end
if ~isempty(v_B)
    F = F - sum(v_B(:,2).'.*K(:,v_B(:,1)+N),2);
    index_delete = [index_delete;v_B(:,1)+N];
    v_index_I(v_B(:,1),:) = [];
    v_index_I = [(1:size(v_index_I,1))'+size(u_index_I,1),v_index_I];
end
if ~isempty(w_B)
    F = F - sum(w_B(:,2).'.*K(:,w_B(:,1)+2*N),2);
    index_delete = [index_delete;w_B(:,1)+2*N];
    w_index_I(w_B(:,1),:) = [];
    w_index_I = [(1:size(w_index_I,1))'+size(u_index_I,1)+size(v_index_I,1),w_index_I];
end
if ~isempty(phi_B)
    F = F - sum(phi_B(:,2).'.*K(:,phi_B(:,1)+3*N),2);
    index_delete = [index_delete;phi_B(:,1)+3*N];
    phi_index_I(phi_B(:,1),:) = [];
    phi_index_I = [(1:size(phi_index_I,1))'+size(u_index_I,1)+size(v_index_I,1)+size(w_index_I,1),...
        phi_index_I];
end
K(index_delete,:)=[];
K(:,index_delete)=[];
F(index_delete) = [];
end

