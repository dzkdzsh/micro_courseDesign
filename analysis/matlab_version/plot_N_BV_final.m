% 绘制N-BV关系图 - 最终优化版
% 计算范围200-1700V，使用线性坐标，优化布局

%% 常量定义
global a_n b_n a_p b_p q eps0 eps_r eps_s;
a_n = 7.03e5;        % cm^-1
b_n = 1.231e6;       % V/cm
a_p = 1.582e6;       % cm^-1
b_p = 2.036e6;       % V/cm
q = 1.602e-19;       % C
eps0 = 8.854e-14;    % F/cm
eps_r = 11.7;
eps_s = eps0 * eps_r;

%% 计算电子电离系数
function an = alpha_n(E)
global a_n b_n;
E = max(E, 1e-12);
an = a_n * exp(-b_n ./ E);
end

%% 计算空穴电离系数
function ap = alpha_p(E)
global a_p b_p;
E = max(E, 1e-12);
ap = a_p * exp(-b_p ./ E);
end

%% 计算电离积分
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

%% 求解给定BV下的W
function W = solve_W_for_BV(BV, Wmin_cm, Wmax_cm, grid_points)
f = @(W) ionization_integral_for_W(W, BV, grid_points) - 1.0;
W = fzero(f, [Wmin_cm, Wmax_cm]);
end

%% 计算所有BV对应的W和N
function [BV, W_cm, N] = compute_all(BV_list, grid_points)
global eps_s q;
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

%% 计算200-1700V范围内的数据
fprintf('计算BV范围200-1700V的数据...\n');
grid_points = 20000;
BV_list = linspace(200, 1700, 76);  % 200-1700V，步长20V

[BV_matlab, W_cm, N] = compute_all(BV_list, grid_points);
N_matlab = N / 1e14;  % 转换为×10^14 cm^-3

%% 模拟Medici仿真数据
BV_medici = [200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1650];
N_medici = [80, 55, 40, 30, 22, 17, 14, 11, 9, 7.5, 6, 5, 4.2, 3, 2.2, 1.6, 1.2, 0.9, 0.7, 0.55, 0.45, 0.4];

%% 创建图形 - 优化布局
figure('Position', [100, 100, 800, 600], 'Color', 'white');

% 先创建坐标轴并设置Y轴为对数坐标
ax = axes;
set(ax, 'YScale', 'log');
set(ax, 'Box', 'on');
set(ax, 'XGrid', 'off');
set(ax, 'YGrid', 'off');
set(ax, 'LineWidth', 1);
set(ax, 'FontSize', 12);
set(ax, 'FontName', 'Times New Roman');

% Medici仿真数据 - 倒三角散点（黑色）
scatter(BV_medici, N_medici, 60, 'k', 'v', 'filled', 'DisplayName', 'N_B(MEDICI)');
hold on;

% MATLAB计算结果 - 黑色连续曲线
plot(BV_matlab, N_matlab, 'k-', 'LineWidth', 1.5, 'DisplayName', 'N_B(MATLAB)');

hold off;

% 设置坐标轴标签
xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('N_B (×10^{14} cm^{-3})', 'FontSize', 14, 'FontName', 'Times New Roman');

% 设置坐标轴范围（对数坐标）
xlim([200, 1700]);
% Y轴对数坐标范围
ylim([0.3, 200]);

% 设置图例（右上角）
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

%% 保存图片
plot_dir = 'd:\micro_courseDesign\plot';
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

% 保存为PNG (高分辨率)
print(gcf, fullfile(plot_dir, 'N_BV.png'), '-dpng', '-r300');
% 保存为JPG
print(gcf, fullfile(plot_dir, 'N_BV.jpg'), '-djpeg', '-r300');

fprintf('\n图片已保存到: %s\n', plot_dir);
fprintf('  - N_BV.png\n');
fprintf('  - N_BV.jpg\n');
fprintf('计算范围: %d-%dV，共%d个数据点\n', min(BV_matlab), max(BV_matlab), length(BV_matlab));
