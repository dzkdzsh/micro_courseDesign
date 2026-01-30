% 实现分段线性变换非均匀拉伸方法

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

%% 分段线性变换非均匀拉伸方法
fprintf('=== 分段线性变换非均匀拉伸方法 ===\n');

% 定义分段点（将BV范围分为3个区间）
break_points = [600, 1200];  % BV值的分段点

% 对Medici数据进行分段
segment1_idx = find(BV_medici <= break_points(1));
segment2_idx = find(BV_medici > break_points(1) & BV_medici <= break_points(2));
segment3_idx = find(BV_medici > break_points(2));

% 对每个分段计算Matlab数据和Medici数据的范围
% 对Matlab数据进行插值，使其与Medici数据在相同的BV点上
lgN_matlab_interp = interp1(BV_matlab, lgN_matlab, BV_medici, 'linear');

% 计算每个分段的范围
seg1_matlab_min = min(lgN_matlab_interp(segment1_idx));
seg1_matlab_max = max(lgN_matlab_interp(segment1_idx));
seg1_medici_min = min(lgN_medici_original(segment1_idx));
seg1_medici_max = max(lgN_medici_original(segment1_idx));

seg2_matlab_min = min(lgN_matlab_interp(segment2_idx));
seg2_matlab_max = max(lgN_matlab_interp(segment2_idx));
seg2_medici_min = min(lgN_medici_original(segment2_idx));
seg2_medici_max = max(lgN_medici_original(segment2_idx));

seg3_matlab_min = min(lgN_matlab_interp(segment3_idx));
seg3_matlab_max = max(lgN_matlab_interp(segment3_idx));
seg3_medici_min = min(lgN_medici_original(segment3_idx));
seg3_medici_max = max(lgN_medici_original(segment3_idx));

% 对每个分段进行线性拉伸
lgN_medici_piecewise = zeros(size(lgN_medici_original));

% 分段1: BV <= 600 V
if ~isempty(segment1_idx)
    seg1_scale = (seg1_matlab_max - seg1_matlab_min) / (seg1_medici_max - seg1_medici_min);
    lgN_medici_piecewise(segment1_idx) = (lgN_medici_original(segment1_idx) - seg1_medici_min) * seg1_scale + seg1_matlab_min;
    fprintf('分段1 (BV <= 600 V): 拉伸系数 = %.4f\n', seg1_scale);
end

% 分段2: 600 < BV <= 1200 V
if ~isempty(segment2_idx)
    seg2_scale = (seg2_matlab_max - seg2_matlab_min) / (seg2_medici_max - seg2_medici_min);
    lgN_medici_piecewise(segment2_idx) = (lgN_medici_original(segment2_idx) - seg2_medici_min) * seg2_scale + seg2_matlab_min;
    fprintf('分段2 (600 < BV <= 1200 V): 拉伸系数 = %.4f\n', seg2_scale);
end

% 分段3: BV > 1200 V
if ~isempty(segment3_idx)
    seg3_scale = (seg3_matlab_max - seg3_matlab_min) / (seg3_medici_max - seg3_medici_min);
    lgN_medici_piecewise(segment3_idx) = (lgN_medici_original(segment3_idx) - seg3_medici_min) * seg3_scale + seg3_matlab_min;
    fprintf('分段3 (BV > 1200 V): 拉伸系数 = %.4f\n', seg3_scale);
end

% 计算均方误差
lgN_medici_interp = interp1(BV_medici, lgN_medici_piecewise, BV_matlab, 'linear');
mse = mean((lgN_matlab - lgN_medici_interp).^2);
fprintf('\n分段线性变换均方误差: %.4f\n', mse);

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

lgN_min = min([lgN_matlab, lgN_medici_piecewise]);
lgN_max = max([lgN_matlab, lgN_medici_piecewise]);
y_ticks = linspace(floor(lgN_min), ceil(lgN_max), 6);
set(ax, 'YTick', y_ticks);
set(ax, 'YLim', [floor(lgN_min) - 0.1, ceil(lgN_max) + 0.1]);
set(ax, 'XLim', [200, 1700]);

% 绘制原始数据和变换后的数据
scatter(BV_medici, lgN_medici_original, 40, 'r', 'o', 'filled', 'DisplayName', '原始Medici数据');
hold on;
scatter(BV_medici, lgN_medici_piecewise, 60, 'k', 'v', 'filled', 'DisplayName', '变换后Medici数据');
plot(BV_matlab, lgN_matlab, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Matlab计算数据');

% 绘制分段线
plot([break_points(1), break_points(1)], [lgN_min - 0.2, lgN_max + 0.2], 'g--', 'LineWidth', 1, 'DisplayName', '分段点');
plot([break_points(2), break_points(2)], [lgN_min - 0.2, lgN_max + 0.2], 'g--', 'LineWidth', 1);

hold off;

xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('lg(N_B) (×10^{14} cm^{-3})', 'FontSize', 14, 'FontName', 'Times New Roman');
title(['分段线性变换优化 (MSE=%.4f)' sprintf('%.4f', mse)], 'FontSize', 16, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

% 保存结果
print(gcf, fullfile(plot_dir, 'lg_piecewise_transform.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'lg_piecewise_transform.jpg'), '-djpeg', '-r300');

fprintf('\n分段线性变换优化结果已保存到: %s\n', plot_dir);
fprintf('  - lg_piecewise_transform.png\n');
fprintf('  - lg_piecewise_transform.jpg\n');

% 绘制分段误差分析
figure('Position', [100, 100, 800, 400], 'Color', 'white');

% 计算原始数据和变换后数据的误差
original_error = lgN_matlab_interp - lgN_medici_original;
transformed_error = lgN_matlab_interp - lgN_medici_piecewise;

subplot(1, 2, 1);
plot(BV_medici, original_error, 'r-', 'LineWidth', 1.5, 'DisplayName', '原始数据误差');
hold on;
plot(BV_medici, transformed_error, 'b-', 'LineWidth', 1.5, 'DisplayName', '变换后数据误差');
plot(BV_medici, zeros(size(BV_medici)), 'k--', 'LineWidth', 1, 'DisplayName', '零误差线');
plot([break_points(1), break_points(1)], [min([original_error, transformed_error])-0.1, max([original_error, transformed_error])+0.1], 'g--', 'LineWidth', 1);
plot([break_points(2), break_points(2)], [min([original_error, transformed_error])-0.1, max([original_error, transformed_error])+0.1], 'g--', 'LineWidth', 1);
hold off;
xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('误差 (Matlab - Medici)', 'FontSize', 14, 'FontName', 'Times New Roman');
title('误差对比', 'FontSize', 16, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');
grid on;

subplot(1, 2, 2);
bar_width = 0.35;
x_pos = 1:3;
original_mse = [mean(original_error(segment1_idx).^2), mean(original_error(segment2_idx).^2), mean(original_error(segment3_idx).^2)];
transformed_mse = [mean(transformed_error(segment1_idx).^2), mean(transformed_error(segment2_idx).^2), mean(transformed_error(segment3_idx).^2)];

bar(x_pos - bar_width/2, original_mse, bar_width, 'r', 'DisplayName', '原始数据MSE');
hold on;
bar(x_pos + bar_width/2, transformed_mse, bar_width, 'b', 'DisplayName', '变换后数据MSE');
hold off;
set(gca, 'XTick', x_pos);
set(gca, 'XTickLabel', {'分段1\n(BV <= 600 V)', '分段2\n(600 < BV <= 1200 V)', '分段3\n(BV > 1200 V)'});
xlabel('分段', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('均方误差', 'FontSize', 14, 'FontName', 'Times New Roman');
title('各分段误差对比', 'FontSize', 16, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');
grid on;

% 调整子图布局
tightlayout;

print(gcf, fullfile(plot_dir, 'piecewise_transform_analysis.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'piecewise_transform_analysis.jpg'), '-djpeg', '-r300');

fprintf('\n分析结果已保存到: %s\n', plot_dir);
fprintf('  - piecewise_transform_analysis.png\n');
fprintf('  - piecewise_transform_analysis.jpg\n');

% 返回最佳变换结果
best_lgN_medici = lgN_medici_piecewise;

fprintf('\n=== 分段线性变换优化完成 ===\n');
fprintf('分段点: %.0f V 和 %.0f V\n', break_points(1), break_points(2));
fprintf('总体均方误差: %.4f\n', mse);
