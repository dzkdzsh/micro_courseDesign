% 实现对数比例调整非均匀拉伸方法

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

% 计算lgN值
lgN_matlab = log10(N_matlab);

%% 模拟Medici数据
BV_medici = [200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1650];
N_medici = [80, 55, 40, 30, 22, 17, 14, 11, 9, 7.5, 6, 5, 4.2, 3, 2.2, 1.6, 1.2, 0.9, 0.7, 0.55, 0.45, 0.4];

% 计算lgN值
lgN_medici_original = log10(N_medici);

% 创建输出目录
plot_dir = 'd:\micro_courseDesign\plot';
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

%% 对数比例调整非均匀拉伸方法
fprintf('=== 对数比例调整非均匀拉伸方法 ===\n');

% 定义比例参数范围
scale_values = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5];
mse_values = zeros(length(scale_values), 1);

% 对每个比例参数进行测试
for i = 1:length(scale_values)
    k = scale_values(i);
    
    % 对数比例调整
    % 方法1: 直接调整对数比例
    lgN_medici_scaled = k * lgN_medici_original;
    
    % 对调整后的数据进行范围调整，使其与Matlab数据范围匹配
    min_scaled = min(lgN_medici_scaled);
    max_scaled = max(lgN_medici_scaled);
    min_matlab = min(lgN_matlab);
    max_matlab = max(lgN_matlab);
    
    lgN_medici_adjusted = (lgN_medici_scaled - min_scaled) * (max_matlab - min_matlab) / (max_scaled - min_scaled) + min_matlab;
    
    % 计算均方误差
    lgN_medici_interp = interp1(BV_medici, lgN_medici_adjusted, BV_matlab, 'linear');
    mse_values(i) = mean((lgN_matlab - lgN_medici_interp).^2);
    
    fprintf('比例参数 k=%.1f: 均方误差 = %.4f\n', k, mse_values(i));
end

% 找到最佳比例参数
[min_mse, best_idx] = min(mse_values);
best_scale = scale_values(best_idx);
fprintf('\n最佳比例参数: k=%.1f, 最小均方误差=%.4f\n', best_scale, min_mse);

% 使用最佳比例参数进行最终变换
k = best_scale;
lgN_medici_scaled = k * lgN_medici_original;

% 范围调整
min_scaled = min(lgN_medici_scaled);
max_scaled = max(lgN_medici_scaled);
min_matlab = min(lgN_matlab);
max_matlab = max(lgN_matlab);
lgN_medici_adjusted = (lgN_medici_scaled - min_scaled) * (max_matlab - min_matlab) / (max_scaled - min_scaled) + min_matlab;

% 绘制结果
figure('Position', [100, 100, 800, 600], 'Color', 'white');
ax = axes;
set(ax, 'YScale', 'linear');
set(ax, 'Box', 'on');
set(ax, 'XGrid', 'off');
set(ax, 'YGrid', 'off');
set(ax, 'LineWidth', 1);
set(ax, 'FontSize', 12);
set(ax, 'FontName', 'Times New Roman');

lgN_min = min([lgN_matlab, lgN_medici_adjusted]);
lgN_max = max([lgN_matlab, lgN_medici_adjusted]);
y_ticks = linspace(floor(lgN_min), ceil(lgN_max), 6);
set(ax, 'YTick', y_ticks);
set(ax, 'YLim', [floor(lgN_min) - 0.1, ceil(lgN_max) + 0.1]);
set(ax, 'XLim', [200, 1700]);

% 绘制原始数据和变换后的数据
scatter(BV_medici, lgN_medici_original, 40, 'r', 'o', 'filled', 'DisplayName', '原始Medici数据');
hold on;
scatter(BV_medici, lgN_medici_adjusted, 60, 'k', 'v', 'filled', 'DisplayName', '变换后Medici数据');
plot(BV_matlab, lgN_matlab, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Matlab计算数据');
hold off;

xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('lg(N_B) (×10^{14} cm^{-3})', 'FontSize', 14, 'FontName', 'Times New Roman');
title(['对数比例调整优化 (比例 k=%.1f, MSE=%.4f)' sprintf('%.1f', best_scale) sprintf('%.4f', min_mse)], 'FontSize', 16, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

% 保存结果
print(gcf, fullfile(plot_dir, 'lg_log_scale_transform.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'lg_log_scale_transform.jpg'), '-djpeg', '-r300');

fprintf('\n对数比例调整优化结果已保存到: %s\n', plot_dir);
fprintf('  - lg_log_scale_transform.png\n');
fprintf('  - lg_log_scale_transform.jpg\n');

% 绘制比例参数与均方误差的关系
figure('Position', [100, 100, 600, 400], 'Color', 'white');
plot(scale_values, mse_values, 'b-', 'LineWidth', 2);
hold on;
plot(best_scale, min_mse, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold off;
xlabel('比例参数 k', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('均方误差', 'FontSize', 14, 'FontName', 'Times New Roman');
title('比例参数对均方误差的影响', 'FontSize', 16, 'FontName', 'Times New Roman');
grid on;
legend('均方误差', '最佳参数', 'Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

print(gcf, fullfile(plot_dir, 'log_scale_transform_analysis.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'log_scale_transform_analysis.jpg'), '-djpeg', '-r300');

fprintf('\n分析结果已保存到: %s\n', plot_dir);
fprintf('  - log_scale_transform_analysis.png\n');
fprintf('  - log_scale_transform_analysis.jpg\n');

% 方法2: 对数偏移调整
fprintf('\n=== 对数偏移调整方法 ===\n');

% 定义偏移参数范围
offset_values = [-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5];
mse_values_offset = zeros(length(offset_values), 1);

% 对每个偏移参数进行测试
for i = 1:length(offset_values)
    offset = offset_values(i);
    
    % 对数偏移调整
    lgN_medici_offset = lgN_medici_original + offset;
    
    % 对调整后的数据进行范围调整，使其与Matlab数据范围匹配
    min_offset = min(lgN_medici_offset);
    max_offset = max(lgN_medici_offset);
    min_matlab = min(lgN_matlab);
    max_matlab = max(lgN_matlab);
    
    lgN_medici_adjusted_offset = (lgN_medici_offset - min_offset) * (max_matlab - min_matlab) / (max_offset - min_offset) + min_matlab;
    
    % 计算均方误差
    lgN_medici_interp_offset = interp1(BV_medici, lgN_medici_adjusted_offset, BV_matlab, 'linear');
    mse_values_offset(i) = mean((lgN_matlab - lgN_medici_interp_offset).^2);
    
    fprintf('偏移参数 offset=%.1f: 均方误差 = %.4f\n', offset, mse_values_offset(i));
end

% 找到最佳偏移参数
[min_mse_offset, best_idx_offset] = min(mse_values_offset);
best_offset = offset_values(best_idx_offset);
fprintf('\n最佳偏移参数: offset=%.1f, 最小均方误差=%.4f\n', best_offset, min_mse_offset);

% 比较两种方法
fprintf('\n=== 方法比较 ===\n');
fprintf('对数比例调整方法: 最佳参数 k=%.1f, MSE=%.4f\n', best_scale, min_mse);
fprintf('对数偏移调整方法: 最佳参数 offset=%.1f, MSE=%.4f\n', best_offset, min_mse_offset);

% 选择最佳方法
if min_mse < min_mse_offset
    fprintf('\n选择对数比例调整方法作为最佳方案\n');
    best_lgN_medici = lgN_medici_adjusted;
    best_method = '对数比例调整';
    best_param = best_scale;
    best_mse = min_mse;
else
    fprintf('\n选择对数偏移调整方法作为最佳方案\n');
    % 使用最佳偏移参数进行最终变换
    lgN_medici_offset = lgN_medici_original + best_offset;
    min_offset = min(lgN_medici_offset);
    max_offset = max(lgN_medici_offset);
    min_matlab = min(lgN_matlab);
    max_matlab = max(lgN_matlab);
    best_lgN_medici = (lgN_medici_offset - min_offset) * (max_matlab - min_matlab) / (max_offset - min_offset) + min_matlab;
    best_method = '对数偏移调整';
    best_param = best_offset;
    best_mse = min_mse_offset;
end

% 绘制最佳方法的结果
figure('Position', [100, 100, 800, 600], 'Color', 'white');
ax = axes;
set(ax, 'YScale', 'linear');
set(ax, 'Box', 'on');
set(ax, 'XGrid', 'off');
set(ax, 'YGrid', 'off');
set(ax, 'LineWidth', 1);
set(ax, 'FontSize', 12);
set(ax, 'FontName', 'Times New Roman');

lgN_min = min([lgN_matlab, best_lgN_medici]);
lgN_max = max([lgN_matlab, best_lgN_medici]);
y_ticks = linspace(floor(lgN_min), ceil(lgN_max), 6);
set(ax, 'YTick', y_ticks);
set(ax, 'YLim', [floor(lgN_min) - 0.1, ceil(lgN_max) + 0.1]);
set(ax, 'XLim', [200, 1700]);

% 绘制原始数据和变换后的数据
scatter(BV_medici, lgN_medici_original, 40, 'r', 'o', 'filled', 'DisplayName', '原始Medici数据');
hold on;
scatter(BV_medici, best_lgN_medici, 60, 'k', 'v', 'filled', 'DisplayName', '变换后Medici数据');
plot(BV_matlab, lgN_matlab, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Matlab计算数据');
hold off;

xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('lg(N_B) (×10^{14} cm^{-3})', 'FontSize', 14, 'FontName', 'Times New Roman');
title([best_method, '优化 (参数=%.1f, MSE=%.4f)' sprintf('%.1f', best_param) sprintf('%.4f', best_mse)], 'FontSize', 16, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

% 保存最佳方法的结果
print(gcf, fullfile(plot_dir, 'lg_best_log_transform.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'lg_best_log_transform.jpg'), '-djpeg', '-r300');

fprintf('\n最佳方法结果已保存到: %s\n', plot_dir);
fprintf('  - lg_best_log_transform.png\n');
fprintf('  - lg_best_log_transform.jpg\n');

fprintf('\n=== 对数比例调整优化完成 ===\n');
fprintf('最佳方法: %s\n', best_method);
fprintf('最佳参数: %.1f\n', best_param);
fprintf('最小均方误差: %.4f\n', best_mse);
