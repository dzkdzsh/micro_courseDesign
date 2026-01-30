% 实现指数函数变换非均匀拉伸方法

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

%% 指数函数变换非均匀拉伸方法
fprintf('=== 指数函数变换非均匀拉伸方法 ===\n');

% 定义指数参数范围
exp_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5];
mse_values = zeros(length(exp_values), 1);

% 对每个指数参数进行测试
for i = 1:length(exp_values)
    k = exp_values(i);
    
    % 指数函数变换
    % 注意：为了保持数据的正性，我们先将数据平移到正数范围，再进行指数变换，最后平移回来
    min_medici = min(lgN_medici_original);
    shift = abs(min_medici) + 0.1;  % 确保所有值都是正数
    
    lgN_medici_shifted = lgN_medici_original + shift;
    lgN_medici_transformed = exp(k * lgN_medici_shifted);
    % 对变换后的数据进行对数变换，使其回到原始尺度附近
    lgN_medici_exp = log(lgN_medici_transformed) / k - shift;
    
    % 对变换后的数据进行范围调整，使其与Matlab数据范围匹配
    min_transformed = min(lgN_medici_exp);
    max_transformed = max(lgN_medici_exp);
    min_matlab = min(lgN_matlab);
    max_matlab = max(lgN_matlab);
    
    lgN_medici_scaled = (lgN_medici_exp - min_transformed) * (max_matlab - min_matlab) / (max_transformed - min_transformed) + min_matlab;
    
    % 计算均方误差
    lgN_medici_interp = interp1(BV_medici, lgN_medici_scaled, BV_matlab, 'linear');
    mse_values(i) = mean((lgN_matlab - lgN_medici_interp).^2);
    
    fprintf('指数参数 k=%.1f: 均方误差 = %.4f\n', k, mse_values(i));
end

% 找到最佳指数参数
[min_mse, best_idx] = min(mse_values);
best_exp = exp_values(best_idx);
fprintf('\n最佳指数参数: k=%.1f, 最小均方误差=%.4f\n', best_exp, min_mse);

% 使用最佳指数参数进行最终变换
k = best_exp;
min_medici = min(lgN_medici_original);
shift = abs(min_medici) + 0.1;
lgN_medici_shifted = lgN_medici_original + shift;
lgN_medici_transformed = exp(k * lgN_medici_shifted);
lgN_medici_exp = log(lgN_medici_transformed) / k - shift;

% 范围调整
min_transformed = min(lgN_medici_exp);
max_transformed = max(lgN_medici_exp);
min_matlab = min(lgN_matlab);
max_matlab = max(lgN_matlab);
lgN_medici_scaled = (lgN_medici_exp - min_transformed) * (max_matlab - min_matlab) / (max_transformed - min_transformed) + min_matlab;

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

lgN_min = min([lgN_matlab, lgN_medici_scaled]);
lgN_max = max([lgN_matlab, lgN_medici_scaled]);
y_ticks = linspace(floor(lgN_min), ceil(lgN_max), 6);
set(ax, 'YTick', y_ticks);
set(ax, 'YLim', [floor(lgN_min) - 0.1, ceil(lgN_max) + 0.1]);
set(ax, 'XLim', [200, 1700]);

% 绘制原始数据和变换后的数据
scatter(BV_medici, lgN_medici_original, 40, 'r', 'o', 'filled', 'DisplayName', '原始Medici数据');
hold on;
scatter(BV_medici, lgN_medici_scaled, 60, 'k', 'v', 'filled', 'DisplayName', '变换后Medici数据');
plot(BV_matlab, lgN_matlab, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Matlab计算数据');
hold off;

xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('lg(N_B) (×10^{14} cm^{-3})', 'FontSize', 14, 'FontName', 'Times New Roman');
title(['指数函数变换优化 (参数 k=%.1f, MSE=%.4f)' sprintf('%.1f', best_exp) sprintf('%.4f', min_mse)], 'FontSize', 16, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

% 保存结果
print(gcf, fullfile(plot_dir, 'lg_exponential_transform.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'lg_exponential_transform.jpg'), '-djpeg', '-r300');

fprintf('\n指数函数变换优化结果已保存到: %s\n', plot_dir);
fprintf('  - lg_exponential_transform.png\n');
fprintf('  - lg_exponential_transform.jpg\n');

% 绘制指数参数与均方误差的关系
figure('Position', [100, 100, 600, 400], 'Color', 'white');
plot(exp_values, mse_values, 'b-', 'LineWidth', 2);
hold on;
plot(best_exp, min_mse, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold off;
xlabel('指数参数 k', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('均方误差', 'FontSize', 14, 'FontName', 'Times New Roman');
title('指数参数对均方误差的影响', 'FontSize', 16, 'FontName', 'Times New Roman');
grid on;
legend('均方误差', '最佳参数', 'Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

print(gcf, fullfile(plot_dir, 'exponential_transform_analysis.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'exponential_transform_analysis.jpg'), '-djpeg', '-r300');

fprintf('\n分析结果已保存到: %s\n', plot_dir);
fprintf('  - exponential_transform_analysis.png\n');
fprintf('  - exponential_transform_analysis.jpg\n');

% 返回最佳变换结果
best_lgN_medici = lgN_medici_scaled;
best_exp_value = best_exp;

fprintf('\n=== 指数函数变换优化完成 ===\n');
fprintf('最佳指数参数: %.1f\n', best_exp_value);
fprintf('最小均方误差: %.4f\n', min_mse);
