% 计算拟合结果的脚本
% 本脚本计算并输出三个拟合式子的具体数值

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

% 主计算和拟合
BV_list = linspace(400, 1600, 41);
[BV, W_cm, N] = compute_all(BV_list, 20000);
W_um = W_cm * 1e4;

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

% 输出拟合结果
disp('========================================');
disp('MATLAB拟合结果 (全区间 400-1600 V)');
disp('========================================');
disp(' ');
disp('1. N与BV的关系:');
disp(['   N ≈ ', num2str(C_N, '%.5g'), ' * BV^(', num2str(p_N(1), '%.6f'), ') cm^-3']);
disp(['   (指数 p = ', num2str(p_N(1), '%.6f'), ')']);
disp(' ');
disp('2. W与BV的关系:');
disp(['   W ≈ ', num2str(C_W, '%.5g'), ' * BV^(', num2str(p_W(1), '%.6f'), ') μm']);
disp(['   (指数 p = ', num2str(p_W(1), '%.6f'), ')']);
disp(' ');
disp('3. W与N的关系:');
disp(['   W ≈ ', num2str(C_NW, '%.5g'), ' * N^(', num2str(p_NW(1), '%.6f'), ') μm']);
disp(['   (指数 p = ', num2str(p_NW(1), '%.6f'), ')']);
disp(' ');

% 不同区间的拟合
ranges = {[400,1600], [400,1200], [400,1000], [500,1600], [600,1600], [400,800]};
disp('不同区间的N-BV拟合结果:');
disp('----------------------------------------');
for i = 1:length(ranges)
    range = ranges{i};
    mask = (BV >= range(1)) & (BV <= range(2));
    if any(mask)
        p_range = polyfit(log_BV(mask), log_N(mask), 1);
        C_range = exp(p_range(2));
        disp(['Range ', num2str(range(1)), '-', num2str(range(2)), ' V:']);
        disp(['  N ≈ ', num2str(C_range, '%.6g'), ' * BV^(', num2str(p_range(1), '%.6f'), ')']);
        disp(['  (指数 = ', num2str(p_range(1), '%.6f'), ')']);
        disp(' ');
    end
end

% 保存数据
data = [BV; W_um; N]';
writematrix(data, 'BV_W_N_fits.csv', 'Delimiter', ',');
disp('数据已保存到 BV_W_N_fits.csv');
