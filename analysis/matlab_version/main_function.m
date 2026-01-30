% 主函数实现
% 本文件实现主计算逻辑，包括求解W、计算N、数据拟合等功能

% 常量定义
global a_n b_n a_p b_p q eps0 eps_r eps_s;
a_n = 7.03e5;        % cm^-1
b_n = 1.231e6;       % V/cm
a_p = 1.582e6;       % cm^-1
b_p = 2.036e6;       % V/cm
q = 1.602e-19;       % C
eps0 = 8.854e-14;    % F/cm
eps_r = 11.7;
eps_s = eps0 * eps_r;

% 计算电子电离系数
function an = alpha_n(E)
global a_n b_n;
E = max(E, 1e-12);  % 避免除以零
an = a_n * exp(-b_n ./ E);
end

% 计算空穴电离系数
function ap = alpha_p(E)
global a_p b_p;
E = max(E, 1e-12);  % 避免除以零
ap = a_p * exp(-b_p ./ E);
end

% 计算电离积分
function I = ionization_integral_for_W(W_cm, BV, grid_points)
% 生成空间网格（确保是列向量）
x = linspace(0, W_cm, grid_points)';
% 计算电场分布
E = 2.0 * BV * (W_cm - x) / (W_cm^2);
% 计算电离系数
an = alpha_n(E);
% 确保an是列向量
an = an(:);
% 计算差值
diff = an - alpha_p(E);
% 确保diff是列向量
diff = diff(:);
% 累积积分
S = cumtrapz(x, diff);
% 确保S是列向量
S = S(:);
% 确保S的长度与x相同
if length(S) < length(x)
    S = [0; S];
end
% 截断或填充S，确保长度与an相同
S = S(1:length(an));
% 计算被积函数
integrand = an .* exp(-S);
% 计算积分
I = trapz(x, integrand);
end

% 求解给定BV下的W
function W = solve_W_for_BV(BV, Wmin_cm, Wmax_cm, grid_points)
% 定义目标函数
f = @(W) ionization_integral_for_W(W, BV, grid_points) - 1.0;
% 使用fzero求解
W = fzero(f, [Wmin_cm, Wmax_cm]);
end

% 计算所有BV对应的W和N
% 输入参数：
%   BV_list - 击穿电压列表 (V)，默认值为400-1600V，步长30V
%   grid_points - 空间网格点数，默认值为20000
% 输出参数：
%   BV - 击穿电压列表
%   W_cm - 耗尽区宽度列表 (cm)
%   N - 掺杂浓度列表 (cm^-3)
function [BV, W_cm, N] = compute_all(varargin)
global eps_s q;
% 设置默认参数
if nargin < 1
    BV_list = linspace(400, 1600, 41);  % 默认BV范围
else
    BV_list = varargin{1};
end

if nargin < 2
    grid_points = 20000;  % 默认网格点数
else
    grid_points = varargin{2};
end

% 参数检查
if ~isvector(BV_list)
    error('BV_list must be a vector');
end

if grid_points <= 0
    error('grid_points must be positive');
end

W_solutions = [];
N_solutions = [];

for i = 1:length(BV_list)
    current_BV = BV_list(i);
    % 使用Fulop近似作为初始猜测
    W_fulop_um = 0.0257 * (current_BV^(7.0/6.0));
    W_guess_cm = W_fulop_um * 1e-4;
    Wmin = max(1e-6, W_guess_cm * 0.2);
    Wmax = max(Wmin*10, W_guess_cm * 4.0);

    try
        W = solve_W_for_BV(current_BV, Wmin, Wmax, grid_points);
    catch
        % 扩大搜索范围
        W = solve_W_for_BV(current_BV, Wmin*0.05, Wmax*20, grid_points);
    end

    % 计算N
    N_val = 2.0 * eps_s * current_BV / (q * W^2);
    W_solutions = [W_solutions, W];
    N_solutions = [N_solutions, N_val];
end

BV = BV_list;
W_cm = W_solutions;
N = N_solutions;
end

% 主函数
% 功能：执行完整的计算流程，包括求解W、计算N、数据拟合和结果输出
% 使用方法：
%   1. 在MATLAB命令窗口中直接运行：main()
%   2. 或调用compute_all函数获取计算结果
function main()
% 设置参数
grid_points = 20000;
BV_list = linspace(400, 1600, 41);

% 计算W和N
[BV, W_cm, N] = compute_all(BV_list, grid_points);

% 转换单位
W_um = W_cm * 1e4;  % 转换为μm

% 拟合N与BV的关系
log_BV = log(BV);
log_N = log(N);
p_N = polyfit(log_BV, log_N, 1);
C_N = exp(p_N(2));

% 拟合W与BV的关系
log_W = log(W_um);
p_W = polyfit(log_BV, log_W, 1);
C_W = exp(p_W(2));

% 拟合W与N的关系
log_N_vals = log(N);
p_NW = polyfit(log_N_vals, log_W, 1);
C_NW = exp(p_NW(2));

% 输出结果
disp('拟合结果 (MATLAB数值拟合):');
disp(['W ≈ ', num2str(C_W, '%.5g'), ' * BV^', num2str(p_W(1), '%.5f'), ' μm']);
disp(['N ≈ ', num2str(C_N, '%.5g'), ' * BV^', num2str(p_N(1), '%.5f'), ' cm^-3']);
disp(['W ≈ ', num2str(C_NW, '%.5g'), ' * N^', num2str(p_NW(1), '%.5f'), ' μm']);

% 保存数据
data = [BV; W_um; N]';
writematrix(data, 'BV_W_N.mat', 'Delimiter', ' ');
writematrix(data, 'BV_W_N.csv', 'Delimiter', ',');

% 不同区间的拟合
ranges = {[400,1600], [400,1200], [400,1000], [500,1600], [600,1600], [400,800]};
disp('\n不同区间的拟合结果:');
for i = 1:length(ranges)
    range = ranges{i};
    mask = (BV >= range(1)) & (BV <= range(2));
    if any(mask)
        p_range = polyfit(log_BV(mask), log_N(mask), 1);
        C_range = exp(p_range(2));
        disp(['range ', num2str(range(1)), '-', num2str(range(2)), ': exponent=', num2str(p_range(1), '%.6f'), ', C=', num2str(C_range, '%.6g')]);
    end
end

% 使用示例
disp('\n使用示例:');
disp('1. 运行完整计算: main()');
disp('2. 计算特定BV范围: [BV, W, N] = compute_all([500, 1000, 1500], 10000);');
disp('3. 使用默认参数: [BV, W, N] = compute_all();');
end