% 计算 R_on,pp(BV) 的脚本
% 本脚本计算并输出导通电阻与击穿电压的关系
%
% 数学模型：
%   R_on,pp = W_pp / (q × μ_n × N_pp)
%
% 其中电子迁移率 μ_n 采用 Caughey-Thomas 模型 [14, 15]：
%   μ_n = 55.24 + 1373.99 / (1 + (N / (1.072 × 10^17))^0.73)
%
% 参考文献：
%   [14] Caughey, D. M., & Thomas, R. E. (1967). Carrier mobilities in silicon
%        empirically related to doping and field. Proceedings of the IEEE,
%        55(12), 2192-2193.
%   [15] Arora, N. D., Hauser, J. R., & Roulston, D. J. (1982). Electron and
%        hole mobilities in silicon as a function of concentration and temperature.
%        IEEE Transactions on Electron Devices, 29(2), 292-295.

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

% Caughey-Thomas 模型参数
global mu_n0 mu_nmax N_ref alpha_CT;
mu_n0 = 55.24;           % cm^2/(V·s)，低场迁移率极限
mu_nmax = 1373.99;       % cm^2/(V·s)，高场迁移率极限
N_ref = 1.072e17;        % cm^-3，参考掺杂浓度
alpha_CT = 0.73;         % 幂指数

% 计算电子电离系数
function an = alpha_n(E)
global a_n b_n;
E = max(E, 1e-12);
an = a_n * exp(-b_n ./ E);
end

% 计算空穴电离系数
function ap = alpha_p(E)
global a_p b_p;
E = max(E, 1e-12);
ap = a_p * exp(-b_p ./ E);
end

% 计算电离积分
function I = ionization_integral_for_W(W_cm, BV, grid_points)
x = linspace(0, W_cm, grid_points)';
E = 2.0 * BV * (W_cm - x) / (W_cm^2);
an = alpha_n(E);
an = an(:);
diff = an - alpha_p(E);
diff = diff(:);
S = cumtrapz(x, diff);
S = S(:);
if length(S) < length(x)
    S = [0; S];
end
S = S(1:length(an));
integrand = an .* exp(-S);
I = trapz(x, integrand);
end

% 求解给定BV下的W
function W = solve_W_for_BV(BV, Wmin_cm, Wmax_cm, grid_points)
f = @(W) ionization_integral_for_W(W, BV, grid_points) - 1.0;
W = fzero(f, [Wmin_cm, Wmax_cm]);
end

% 计算所有BV对应的W和N
function [BV, W_cm, N] = compute_all(BV_list, grid_points)
global eps_s q;
W_solutions = [];
N_solutions = [];

for i = 1:length(BV_list)
    current_BV = BV_list(i);
    W_fulop_um = 0.0257 * (current_BV^(7.0/6.0));
    W_guess_cm = W_fulop_um * 1e-4;
    Wmin = max(1e-6, W_guess_cm * 0.2);
    Wmax = max(Wmin*10, W_guess_cm * 4.0);

    try
        W = solve_W_for_BV(current_BV, Wmin, Wmax, grid_points);
    catch
        W = solve_W_for_BV(current_BV, Wmin*0.05, Wmax*20, grid_points);
    end

    N_val = 2.0 * eps_s * current_BV / (q * W^2);
    W_solutions = [W_solutions, W];
    N_solutions = [N_solutions, N_val];
end

BV = BV_list;
W_cm = W_solutions;
N = N_solutions;
end

% 计算电子迁移率（Caughey-Thomas 模型）
function mu_n = electron_mobility_ct(N)
% Caughey-Thomas 模型计算电子迁移率
% 输入：N - 掺杂浓度 (cm^-3)
% 输出：mu_n - 电子迁移率 (cm^2/(V·s))
%
% 模型公式：
%   μ_n = μ_n0 + μ_nmax / (1 + (N / N_ref)^alpha)
%
% 其中：
%   μ_n0 = 55.24 cm^2/(V·s)      - 低场迁移率极限
%   μ_nmax = 1373.99 cm^2/(V·s)  - 高场迁移率极限
%   N_ref = 1.072 × 10^17 cm^-3  - 参考掺杂浓度
%   alpha = 0.73                 - 幂指数

global mu_n0 mu_nmax N_ref alpha_CT;

% 确保N为正数
N = max(N, 1e-10);

% Caughey-Thomas 公式
mu_n = mu_n0 + mu_nmax ./ (1 + (N ./ N_ref).^alpha_CT);
end

% 计算导通电阻 R_on,pp
function R_on = calculate_Ron_pp(W_cm, N)
% 计算导通电阻 R_on,pp
% 输入：
%   W_cm - 耗尽区宽度 (cm)
%   N    - 掺杂浓度 (cm^-3)
% 输出：
%   R_on - 导通电阻 (Ω·cm^2)
%
% 公式：R_on,pp = W_pp / (q × μ_n × N_pp)
%
% 注意：W_cm单位为cm，结果单位为 Ω·cm^2

global q;

% 计算电子迁移率
mu_n = electron_mobility_ct(N);

% 计算导通电阻
% R_on = W / (q * μ_n * N)
% 单位：cm / (C * cm^2/(V·s) * cm^-3) = cm / (C/(V·s·cm)) = cm * V·s·cm / C = V·s·cm^2 / C = Ω·cm^2
R_on = W_cm ./ (q .* mu_n .* N);
end

% 主计算
BV_list = linspace(400, 1600, 41);
[BV, W_cm, N] = compute_all(BV_list, 20000);
W_um = W_cm * 1e4;

% 计算电子迁移率
mu_n = electron_mobility_ct(N);

% 计算导通电阻
R_on_pp = calculate_Ron_pp(W_cm, N);

% 拟合 R_on,pp 与 BV 的关系
log_BV = log(BV);
log_Ron = log(R_on_pp);
p_Ron = polyfit(log_BV, log_Ron, 1);
C_Ron = exp(p_Ron(2));

% 输出结果
disp('========================================');
disp('R_on,pp(BV) 计算结果');
disp('========================================');
disp(' ');
disp('数学模型：');
disp('  R_on,pp = W_pp / (q × μ_n × N_pp)');
disp(' ');
disp('电子迁移率（Caughey-Thomas 模型）[14, 15]：');
disp('  μ_n = 55.24 + 1373.99 / (1 + (N / (1.072×10^17))^0.73)');
disp(' ');
disp('----------------------------------------');
disp('拟合结果（全区间 400-1600 V）：');
disp(['  R_on,pp ≈ ', num2str(C_Ron, '%.5g'), ' × BV^(', num2str(p_Ron(1), '%.6f'), ') Ω·cm^2']);
disp(['  (指数 p = ', num2str(p_Ron(1), '%.6f'), ')']);
disp(' ');

% 显示部分数据点
disp('部分计算数据（BV, W, N, μ_n, R_on,pp）：');
disp('----------------------------------------');
disp('BV(V)    W(μm)    N(cm^-3)      μ_n(cm^2/Vs)  R_on,pp(Ω·cm^2)');
disp('----------------------------------------------------------------');
indices = [1, 11, 21, 31, 41];  % 显示 400, 700, 1000, 1300, 1600 V
for i = 1:length(indices)
    idx = indices(i);
    fprintf('%4d     %6.3f   %6.3e    %8.2f      %8.4f\n', ...
        BV(idx), W_um(idx), N(idx), mu_n(idx), R_on_pp(idx));
end
disp(' ');

% 不同区间的拟合
ranges = {[400,1600], [400,1200], [400,1000], [500,1600], [600,1600], [400,800]};
disp('不同区间的 R_on,pp-BV 拟合结果：');
disp('----------------------------------------');
for i = 1:length(ranges)
    range = ranges{i};
    mask = (BV >= range(1)) & (BV <= range(2));
    if any(mask)
        p_range = polyfit(log_BV(mask), log_Ron(mask), 1);
        C_range = exp(p_range(2));
        disp(['Range ', num2str(range(1)), '-', num2str(range(2)), ' V:']);
        disp(['  R_on,pp ≈ ', num2str(C_range, '%.6g'), ' × BV^(', num2str(p_range(1), '%.6f'), ')']);
        disp(['  (指数 = ', num2str(p_range(1), '%.6f'), ')']);
        disp(' ');
    end
end

% 保存数据
data = [BV; W_um; N; mu_n; R_on_pp]';
writematrix(data, 'Ron_pp_results.csv', 'Delimiter', ',');
disp('数据已保存到 Ron_pp_results.csv');
disp('列：BV(V), W(μm), N(cm^-3), μ_n(cm^2/Vs), R_on,pp(Ω·cm^2)');
