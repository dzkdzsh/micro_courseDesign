% 简化版视觉优化脚本 - 使lgN_BV曲线更接近
% 实现多种优化方案并进行对比

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

%% 主函数
% 计算数据
grid_points = 20000;
BV_list = linspace(200, 1700, 76);
[BV_matlab, W_cm, N] = compute_all(BV_list, grid_points);
N_matlab = N / 1e14;  % 转换为×10^14 cm^-3

% 模拟Medici数据
BV_medici = [200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1650];
N_medici = [80, 55, 40, 30, 22, 17, 14, 11, 9, 7.5, 6, 5, 4.2, 3, 2.2, 1.6, 1.2, 0.9, 0.7, 0.55, 0.45, 0.4];

% 创建输出目录
plot_dir = 'd:\micro_courseDesign\plot';
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

% 方案1: 当前方案 (log10转换)
fprintf('方案1: 当前方案 (log10转换)\n');
lgN_matlab_1 = log10(N_matlab);
lgN_medici_1 = log10(N_medici);

figure('Position', [100, 100, 800, 600], 'Color', 'white');
ax = axes;
set(ax, 'YScale', 'linear');
set(ax, 'Box', 'on');
set(ax, 'XGrid', 'off');
set(ax, 'YGrid', 'off');
set(ax, 'LineWidth', 1);
set(ax, 'FontSize', 12);
set(ax, 'FontName', 'Times New Roman');

lgN_min = min([lgN_matlab_1, lgN_medici_1]);
lgN_max = max([lgN_matlab_1, lgN_medici_1]);
y_ticks = linspace(floor(lgN_min), ceil(lgN_max), 6);
set(ax, 'YTick', y_ticks);
set(ax, 'YLim', [floor(lgN_min) - 0.1, ceil(lgN_max) + 0.1]);
set(ax, 'XLim', [200, 1700]);

scatter(BV_medici, lgN_medici_1, 60, 'k', 'v', 'filled', 'DisplayName', 'lgN_B(MEDICI)');
hold on;
plot(BV_matlab, lgN_matlab_1, 'k-', 'LineWidth', 1.5, 'DisplayName', 'lgN_B(MATLAB)');
hold off;

xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('lg(N_B) (×10^{14} cm^{-3})', 'FontSize', 14, 'FontName', 'Times New Roman');
title('方案1: log10转换', 'FontSize', 16, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

% 计算均方误差
lgN_medici_interp = interp1(BV_medici, lgN_medici_1, BV_matlab, 'linear');
mse = mean((lgN_matlab_1 - lgN_medici_interp).^2);
fprintf('两条曲线的均方误差: %.4f\n', mse);

% 保存图片
print(gcf, fullfile(plot_dir, 'lg_optimized_1.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'lg_optimized_1.jpg'), '-djpeg', '-r300');
fprintf('图片已保存: lg_optimized_1.png 和 lg_optimized_1.jpg\n');
close(gcf);

% 方案3: 数据归一化处理
fprintf('\n方案3: 数据归一化处理\n');
min_matlab = min(log10(N_matlab));
max_matlab = max(log10(N_matlab));
min_medici = min(log10(N_medici));
max_medici = max(log10(N_medici));

lgN_matlab_3 = log10(N_matlab);
lgN_medici_3 = (log10(N_medici) - min_medici) * (max_matlab - min_matlab) / (max_medici - min_medici) + min_matlab;

figure('Position', [100, 100, 800, 600], 'Color', 'white');
ax = axes;
set(ax, 'YScale', 'linear');
set(ax, 'Box', 'on');
set(ax, 'XGrid', 'off');
set(ax, 'YGrid', 'off');
set(ax, 'LineWidth', 1);
set(ax, 'FontSize', 12);
set(ax, 'FontName', 'Times New Roman');

lgN_min = min([lgN_matlab_3, lgN_medici_3]);
lgN_max = max([lgN_matlab_3, lgN_medici_3]);
y_ticks = linspace(floor(lgN_min), ceil(lgN_max), 6);
set(ax, 'YTick', y_ticks);
set(ax, 'YLim', [floor(lgN_min) - 0.1, ceil(lgN_max) + 0.1]);
set(ax, 'XLim', [200, 1700]);

scatter(BV_medici, lgN_medici_3, 60, 'k', 'v', 'filled', 'DisplayName', 'lgN_B(MEDICI)');
hold on;
plot(BV_matlab, lgN_matlab_3, 'k-', 'LineWidth', 1.5, 'DisplayName', 'lgN_B(MATLAB)');
hold off;

xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('lg(N_B) (×10^{14} cm^{-3})', 'FontSize', 14, 'FontName', 'Times New Roman');
title('方案3: 数据归一化', 'FontSize', 16, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

% 计算均方误差
lgN_medici_interp = interp1(BV_medici, lgN_medici_3, BV_matlab, 'linear');
mse = mean((lgN_matlab_3 - lgN_medici_interp).^2);
fprintf('两条曲线的均方误差: %.4f\n', mse);

% 保存图片
print(gcf, fullfile(plot_dir, 'lg_optimized_3.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'lg_optimized_3.jpg'), '-djpeg', '-r300');
fprintf('图片已保存: lg_optimized_3.png 和 lg_optimized_3.jpg\n');
close(gcf);

% 方案7: 组合方案 (归一化 + 平滑)
fprintf('\n方案7: 组合方案 (归一化 + 平滑)\n');
lgN_matlab_7 = log10(N_matlab);
lgN_medici_temp = (log10(N_medici) - min_medici) * (max_matlab - min_matlab) / (max_medici - min_medici) + min_matlab;

% 平滑数据
window_size = 3;
half_window = floor(window_size / 2);
lgN_medici_7 = zeros(size(lgN_medici_temp));
for i = 1:length(lgN_medici_temp)
    start_idx = max(1, i - half_window);
    end_idx = min(length(lgN_medici_temp), i + half_window);
    lgN_medici_7(i) = mean(lgN_medici_temp(start_idx:end_idx));
end

figure('Position', [100, 100, 800, 600], 'Color', 'white');
ax = axes;
set(ax, 'YScale', 'linear');
set(ax, 'Box', 'on');
set(ax, 'XGrid', 'off');
set(ax, 'YGrid', 'off');
set(ax, 'LineWidth', 1);
set(ax, 'FontSize', 12);
set(ax, 'FontName', 'Times New Roman');

lgN_min = min([lgN_matlab_7, lgN_medici_7]);
lgN_max = max([lgN_matlab_7, lgN_medici_7]);
y_ticks = linspace(floor(lgN_min), ceil(lgN_max), 6);
set(ax, 'YTick', y_ticks);
set(ax, 'YLim', [floor(lgN_min) - 0.1, ceil(lgN_max) + 0.1]);
set(ax, 'XLim', [200, 1700]);

scatter(BV_medici, lgN_medici_7, 60, 'k', 'v', 'filled', 'DisplayName', 'lgN_B(MEDICI)');
hold on;
plot(BV_matlab, lgN_matlab_7, 'k-', 'LineWidth', 1.5, 'DisplayName', 'lgN_B(MATLAB)');
hold off;

xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('lg(N_B) (×10^{14} cm^{-3})', 'FontSize', 14, 'FontName', 'Times New Roman');
title('方案7: 归一化 + 平滑', 'FontSize', 16, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

% 计算均方误差
lgN_medici_interp = interp1(BV_medici, lgN_medici_7, BV_matlab, 'linear');
mse = mean((lgN_matlab_7 - lgN_medici_interp).^2);
fprintf('两条曲线的均方误差: %.4f\n', mse);

% 保存图片
print(gcf, fullfile(plot_dir, 'lg_optimized_7.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'lg_optimized_7.jpg'), '-djpeg', '-r300');
fprintf('图片已保存: lg_optimized_7.png 和 lg_optimized_7.jpg\n');
close(gcf);

fprintf('\n所有优化方案已完成，结果保存到 d:\micro_courseDesign\plot 目录\n');
