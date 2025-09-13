% 有精确解的特征值问题算例
clear
format long
% clc
% currentPath = 'D:\Code\M\Mortar_FEM_Wavelet';
currentPath = 'D:\code\Mortar_FEM_Wavelet3D';
addpath(genpath(currentPath));
% 使用的数值方法
% method_type="power";
% method_type="J-D";
method_type="eigs";
% 无量纲化？
dimensionless=true;
% dimensionless=false;
% 材料参数
% ========== 压电材料参数 (PZT-5A 量级) ==========
% 所有参数使用 Voigt 表示法

% 密度 [kg/m³]
rho = 7500;

% 弹性矩阵 (6×6, 单位: Pa)
% Voigt 顺序: 11, 22, 33, 23, 13, 12
c_LN = 1e9 * [ % 乘以 10^9 转换为 Pa
    120,  75,  75,   0,   0,   0;
    75, 120,  75,   0,   0,   0;
    75,  75, 120,   0,   0,   0;
    0,   0,   0, 22.5,  0,   0;
    0,   0,   0,   0, 22.5,  0;
    0,   0,   0,   0,   0, 22.5
    ];

% 压电矩阵 (3×6, 单位: C/m²)
% 行: 电场方向 (1=x, 2=y, 3=z)
% 列: Voigt 力学索引
e_LN = [
    0,   0,   0,   0, 12.3,  0;
    0,   0,   0, 12.3,   0,  0;
    -5.4,-5.4,15.8,   0,   0,  0
    ];

% 介电矩阵 (3×3, 单位: F/m)
epcl_LN = 1e-9 * [ % 乘以 10^{-9}
    8.85, 0,    0;
    0,    8.85, 0;
    0,    0,    8.85
    ];

% ========== 最小特征值精确解参数 ==========
% 特征频率 [rad/s]
omega_exact = pi * sqrt(c_LN(5,5) / rho); % ≈ 5442 rad/s

% 位移场函数 (x,y,z 为坐标向量)
u_x = @(x,y,z) sin(pi*z);
u_y = @(x,y,z) 0;
u_z = @(x,y,z) 0;

% 电势场函数
phi = @(x,y,z) 0;

% 几何参数
ori = [0,0,0];
a = 1; b = 1; c = 1;
kappa_bar_sq=1;
% 无量纲化
if dimensionless==true
    [c_LN, e_LN, epcl_LN, L_bar, L0, vs, Phi0, kappa_bar_sq] = ...
        piezo_dimensionless(c_LN, e_LN, epcl_LN, rho, a, b, c); rho=1;
end
% 离散
% type="quadratic";
type="linear";
N = 32;
Nx=N+1;Ny=N+1;Nz=N+1;
[K,M,~,Dof_Index] =...
    AssemblePiezMatFEM(c_LN,e_LN,epcl_LN,rho,kappa_bar_sq,...
    ori,a,b,c,Nx,Ny,Nz,...
    [],[],[],type);
% 边界条件
bound_vanish=find(Dof_Index(:,4)==ori(3)|Dof_Index(:,4)==ori(3)+c);
bound_xL = find(Dof_Index(:,2)==ori(1));
bound_xR = find(Dof_Index(:,2)==ori(1)+a);
bound_yL = find(Dof_Index(:,3)==ori(2));
bound_yR = find(Dof_Index(:,3)==ori(2)+b);

K(bound_xL,:)=K(bound_xL,:)+K(bound_xR,:);
K(bound_yL,:)=K(bound_yL,:)+K(bound_yR,:);
K(:,bound_xL)=K(:,bound_xL)+K(:,bound_xR);
K(:,bound_yL)=K(:,bound_yL)+K(:,bound_yR);

M(bound_xL,:)=M(bound_xL,:)+M(bound_xR,:);
M(bound_yL,:)=M(bound_yL,:)+M(bound_yR,:);
M(:,bound_xL)=M(:,bound_xL)+M(:,bound_xR);
M(:,bound_yL)=M(:,bound_yL)+M(:,bound_yR);

bound_vanish=[bound_vanish;bound_xR;bound_yR];

K(bound_vanish,:)=[];
K(:,bound_vanish)=[];
M(bound_vanish,:)=[];
M(:,bound_vanish)=[];
Dof_Index(bound_vanish,:)=[];

%% 计算特征值
tic
if method_type == "power"
    [lambda_saddle_eigs, V_saddle_eigs, res, dist] = inverse_iteration_GPT(K, M, 0, 1e-10, 1000);
elseif method_type =="J-D"
    [lambda_saddle_eigs, V_saddle_eigs, res, dist]=JD_iteration_Gen(K,[],M,1E-10, 200);
elseif method_type =="eigs"
    [V_saddle_eigs,lambda_saddle_eigs]=eigs(K,M,1,'sm');
end
toc
if dimensionless == true
    [lambda_saddle_eigs, ~, ~] = piezo_restore_dimension(lambda_saddle_eigs, V_saddle_eigs, L0, vs, Phi0);
end
err = abs(lambda_saddle_eigs-omega_exact^2)/omega_exact^2;
fprintf('relative error between numerical sol and exact sol = %.2e\n',err);


%%%%%%%%%%%%%%%%%%%%%%%%% 数值求解器 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 反幂法 %%%%%%%%%%%%%
function [lambda, v, res, dist] = inverse_iteration_GPT(A, B, sigma, tol, maxit)
% 求解广义特征值问题 A u = lambda B u
% res：每步的残差|Au^k-lambda^kBu^k|
% dist：相邻两步迭代的变化：|lambda^{k+1}-lambda^k|/|lambda^k|
% 使用反幂(shift-invert)法, 收敛到最靠近 sigma 的特征值
%
% 输入:
%   sigma  - 移位参数
%   tol    - 收敛容忍度 (如 1e-10)
%   maxit  - 最大迭代次数
%
% 输出:
%   lambda - 近似特征值
%   v      - 对应特征向量 (B-单位化)

if nargin < 5, maxit = 200; end
if nargin < 4, tol = 1e-10; end
n = size(A,1);

% 初始向量
% v = randn(n,1);
v=ones(n,1);
v = v / sqrt(v'*v);  % B-归一化
% v = v / sqrt(v'*v);
lambda_old = (v' * A * v) / (v' * B * v);
res=norm(A*v-lambda_old*B*v);
dist=[];
if res(end) < tol
    fprintf('Converged in 0 iterations,\n residual = %.2e\n, distance between two adjoint step is %.2e\n', res);
    return;
end
% 因式分解 A - sigma B
[L,U,P,Q] = lu(A - sigma*B);
for k = 1:maxit
    % 1) 解线性方程 (A - sigma B) w = B v
    rhs = B*v;
    w = Q * (U \ (L \ (P * rhs)));   % 使用LU分解求解

    % 2) B-单位化
    v=w/sqrt(w'*w);
    % 3) Rayleigh 商 (近似特征值)
    lambda = (v' * A * v) / (v' * B * v);

    % 4) 收敛判据
    res = [res;norm(A*v - lambda*B*v)];
    dist = [dist;(lambda - lambda_old)/abs(lambda_old)];
    if res(end) < tol || abs(dist(end)) < tol
        fprintf('Converged in %d iterations,\n residual = %.2e\n, distance between two adjoint step is %.2e\n', k, res(end),dist(end));
        return;
    end
    lambda_old = lambda;
end

warning('未收敛，最大迭代数 = %d, 最后残差 = %.2e', maxit, res);
end

%%%%%%%%Jacobi-Davidson 方法%%%%%%%%%%
function [lambda,vec,res,dist]=JD_iteration_Gen(A,PA,B,tol, maxit)
if ~exist('PA','var')||isempty(PA), PA=A; end
n = size(A,1);
% v = randn(n,1);
v = ones(n,1);
v = v/norm(v);
w=A*v;
w_tilde=B*v;
vAv=v'*w; % Rayleigh商
vBv=v'*w_tilde;
u=v; theta = vAv/vBv;
res=w-theta*w_tilde;
dist=[];
if res<tol
    lambda=theta;
    vec=u;
    fprintf('Converged in 0 iterations, residual = %.2e\n', res);
    return;
end
for k=1:maxit
    lambda_old=theta;
    M = PA-theta*B;
    Mu=B*u; Mu = M\Mu;
    epcl=(u'*u)/(u'*Mu);
    t=-u+epcl*Mu;
    t = modifiedGramSchmidt(t, v);
    t = t/norm(t);
    v=[v,t];
    w=[w,A*t];
    vAv=v'*A*v;
    vBv=v'*B*v;
    [s,theta]=eig(vAv,vBv,'vector');
    [~, idx] = min((theta)); % 以最小模特征值为例
    theta = theta(idx);        % 当前近似特征值
    s = s(:, idx);
    s = s/norm(s);
    u=v*s;
    uhat=w*s;
    r=uhat-theta*B*u;
    res=[res;norm(r)];
    dist=[dist;(theta-lambda_old)/abs(lambda_old)];
    if res(end)<tol || abs(dist(end))<tol
        lambda=theta;
        vec=u;
        fprintf('Converged in %d iterations, residual = %.2e, distance between two adjoint step is %.2e\n', k, res(end),dist(end));
        return;
    end
end
fprintf('Reached the maximum number of %d iterations, residual = %.2e, distance between two adjoint step is %.2e\n', maxit, res(end),dist(end));
lambda=theta;
vec=u;
    function t_orth = modifiedGramSchmidt(t, V)
        % 函数：modifiedGramSchmidt
        % 输入：
        %   t - 待正交化的向量（列向量）
        %   V - 标准正交基组成的矩阵（每列是一个基向量）
        % 输出：
        %   t_orth - 正交化后的向量，满足与V的每一列正交

        % 检查输入维度
        [m, n] = size(V);
        if size(t, 1) ~= m
            error('向量t的维度和矩阵V的行数不一致');
        end
        if size(t, 2) ~= 1
            error('向量t必须是列向量');
        end

        % 初始化正交化结果
        t_orth = t;

        % Modified Gram-Schmidt过程
        for j = 1:n
            vj = V(:, j);          % 取出第j个标准正交基
            proj_coeff = vj' * t_orth; % 计算投影系数
            t_orth = t_orth - proj_coeff * vj; % 减去投影分量
        end
    end
end