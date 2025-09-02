function [c_bar, e_bar, epsilon_bar, L_bar, L0, vs, Phi0, kappa_bar_sq] = ...
         piezo_dimensionless(c_E, e, epsilon_S, rho, Lx, Ly, Lz)
% 压电问题无量纲化程序
% 输入:
%   c_E: 6x6弹性常数矩阵 (Pa)
%   e: 3x6压电应力常数矩阵 (C/m²)
%   epsilon_S: 3x3介电常数矩阵 (F/m)
%   rho: 密度 (kg/m³)
%   Lx, Ly, Lz: 长方体尺寸 (m)
% 输出:
%   c_bar: 无量纲弹性常数矩阵
%   e_bar: 无量纲压电常数矩阵
%   epsilon_bar: 无量纲介电常数矩阵
%   L_bar: 无量纲尺寸 [Lx_bar, Ly_bar, Lz_bar]
%   L0: 特征长度 (m)
%   vs: 特征声速 (m/s)
%   Phi0: 特征电势 (V)
%   kappa_bar_sq: 无量纲机电耦合系数

% 1. 计算特征长度 (取最小尺寸)
L0 = min([Lx, Ly, Lz]);

% 2. 计算特征声速 (基于最大弹性模量)
% 提取弹性矩阵的主对角线元素(简化处理)
c_diag = diag(c_E);
c_max = max(c_diag(1:3)); % 取xx,yy,zz方向的最大值
vs = sqrt(c_max / rho);

% 3. 计算特征压电常数 (取矩阵中的最大绝对值)
e0 = max(abs(e(:)));

% 4. 计算特征介电常数 (取矩阵中的最大绝对值)
epsilon0 = max(abs(diag(epsilon_S))); % 通常介电常数矩阵对角占优

% 5. 计算特征电势
Phi0 = (e0 * L0) / epsilon0;

% 6. 计算无量纲机电耦合系数
kappa_bar_sq = e0^2 / (epsilon0 * rho * vs^2);

% 7. 无量纲化材料参数
scale_stress = rho * vs^2; % 应力缩放因子
c_bar = c_E / scale_stress;
e_bar = e / e0;
epsilon_bar = epsilon_S / epsilon0;

% 8. 无量纲化几何尺寸
L_bar = [Lx, Ly, Lz] / L0;

% 显示无量纲化结果
fprintf('无量纲化完成:\n');
fprintf('特征长度 L0 = %.4e m\n', L0);
fprintf('特征声速 vs = %.4e m/s\n', vs);
fprintf('特征电势 Phi0 = %.4e V\n', Phi0);
fprintf('无量纲机电耦合系数 kappa_bar^2 = %.4f\n', kappa_bar_sq);
fprintf('无量纲尺寸: Lx_bar = %.2f, Ly_bar = %.2f, Lz_bar = %.2f\n', L_bar);
end
