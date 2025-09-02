function [omega, U, phi] = piezo_restore_dimension(omega_bar, U_bar, L0, vs, Phi0)
% 压电问题量纲还原程序
% 输入:
%   omega_bar: 无量纲特征频率 (标量或向量)
%   U_bar: 无量纲特征函数矩阵 [N x M]
%         每列是一个特征向量，排列为 [u1_bar, u2_bar, u3_bar, phi_bar]^T
%   L0: 特征长度 (m)
%   vs: 特征声速 (m/s)
%   Phi0: 特征电势 (V)
% 输出:
%   omega: 有量纲特征频率 (rad/s)
%   U: 有量纲位移特征函数 [N/4 x 3 x M]
%   phi: 有量纲电势特征函数 [N/4 x M]

% 1. 还原特征频率
omega = omega_bar * (vs / L0)^2;

% 2. 确定特征函数数量
num_modes = size(U_bar, 2);
num_dofs = size(U_bar, 1);

% 3. 检查自由度数量 (应为4的倍数)
if mod(num_dofs, 4) ~= 0
    error('自由度数量必须是4的倍数 (每个节点有ux,uy,uz,phi四个自由度)');
end

num_nodes = num_dofs / 4;

% 4. 初始化输出数组
U = zeros(num_nodes, 3, num_modes);
phi = zeros(num_nodes, num_modes);

% 5. 还原每个特征函数
for m = 1:num_modes
    % 提取当前模态的无量纲特征向量
    u_bar = U_bar(:, m);
    
    % 还原位移分量 (乘以L0)
    U(:, 1, m) = u_bar(1:4:end) * L0; % x-位移
    U(:, 2, m) = u_bar(2:4:end) * L0; % y-位移
    U(:, 3, m) = u_bar(3:4:end) * L0; % z-位移
    
    % 还原电势分量 (乘以Phi0)
    phi(:, m) = u_bar(4:4:end) * Phi0;
end

fprintf('量纲还原完成:\n');
%fprintf('特征频率范围: %.2f ~ %.2f rad/s\n', min(omega), max(omega));
end