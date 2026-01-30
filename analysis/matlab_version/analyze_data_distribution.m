% 分析Matlab计算数据和Medici数据的分布特征

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

%% 主分析函数
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
lgN_medici = log10(N_medici);

%% 数据分布特征分析
fprintf('=== 数据分布特征分析 ===\n');

% Matlab数据特征
fprintf('\nMatlab计算数据特征:\n');
fprintf('数据点数量: %d\n', length(lgN_matlab));
fprintf('lgN范围: %.4f 到 %.4f\n', min(lgN_matlab), max(lgN_matlab));
fprintf('lgN平均值: %.4f\n', mean(lgN_matlab));
fprintf('lgN标准差: %.4f\n', std(lgN_matlab));

% Medici数据特征
fprintf('\nMedici数据特征:\n');
fprintf('数据点数量: %d\n', length(lgN_medici));
fprintf('lgN范围: %.4f 到 %.4f\n', min(lgN_medici), max(lgN_medici));
fprintf('lgN平均值: %.4f\n', mean(lgN_medici));
fprintf('lgN标准差: %.4f\n', std(lgN_medici));

% 数据差异分析
fprintf('\n数据差异分析:\n');
% 对Medici数据进行插值，使其与Matlab数据在相同的BV点上
lgN_medici_interp = interp1(BV_medici, lgN_medici, BV_matlab, 'linear');
difference = lgN_matlab - lgN_medici_interp;
fprintf('平均差异: %.4f\n', mean(difference));
fprintf('最大差异: %.4f\n', max(abs(difference)));
fprintf('均方误差: %.4f\n', mean(difference.^2));

% 差异随BV的变化
fprintf('\n差异随BV的变化趋势:\n');
BVT_Bins = linspace(200, 1700, 6);
for i = 1:length(BVT_Bins)-1
    mask = (BV_matlab >= BVT_Bins(i)) & (BV_matlab < BVT_Bins(i+1));
    if any(mask)
        bin_diff = mean(difference(mask));
        fprintf('BV %.0f-%.0f V: 平均差异 = %.4f\n', BVT_Bins(i), BVT_Bins(i+1), bin_diff);
    end
end

% 数据分布可视化
figure('Position', [100, 100, 1000, 600], 'Color', 'white');

% 子图1: 原始数据对比
subplot(2, 2, 1);
scatter(BV_medici, lgN_medici, 60, 'k', 'v', 'filled', 'DisplayName', 'lgN_B(MEDICI)');
hold on;
plot(BV_matlab, lgN_matlab, 'k-', 'LineWidth', 1.5, 'DisplayName', 'lgN_B(MATLAB)');
hold off;
xlabel('BV_{pp} (V)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('lg(N_B) (×10^{14} cm^{-3})', 'FontSize', 12, 'FontName', 'Times New Roman');
title('原始数据对比', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman');
grid on;

% 子图2: 数据差异
subplot(2, 2, 2);
plot(BV_matlab, difference, 'r-', 'LineWidth', 1.2, 'DisplayName', '差异 (Matlab - Medici)');
xlabel('BV_{pp} (V)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('lg(N_B) 差异', 'FontSize', 12, 'FontName', 'Times New Roman');
title('数据差异随BV的变化', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman');
grid on;

% 子图3: 数据分布直方图
subplot(2, 2, 3);
histogram(lgN_matlab, 20, 'FaceColor', 'blue', 'EdgeColor', 'black', 'DisplayName', 'Matlab数据');
hold on;
histogram(lgN_medici, 15, 'FaceColor', 'red', 'EdgeColor', 'black', 'DisplayName', 'Medici数据');
hold off;
xlabel('lg(N_B) (×10^{14} cm^{-3})', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('频次', 'FontSize', 12, 'FontName', 'Times New Roman');
title('数据分布直方图', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman');
grid on;

% 子图4: 对数-对数关系
subplot(2, 2, 4);
plot(lgN_medici, lgN_medici, 'k--', 'LineWidth', 1, 'DisplayName', '理想匹配');
hold on;
% 对Medici数据进行插值，使其与Matlab数据点数量相同
lgN_medici_interp_full = interp1(BV_medici, lgN_medici, BV_matlab, 'linear');
scatter(lgN_medici_interp_full, lgN_matlab, 50, 'g', 'o', 'filled', 'DisplayName', '实际数据点');
hold off;
xlabel('lg(N_B) (Medici)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('lg(N_B) (Matlab)', 'FontSize', 12, 'FontName', 'Times New Roman');
title('数据对应关系', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman');
grid on;

% 调整子图布局
tightlayout;

% 保存分析结果
plot_dir = 'd:\micro_courseDesign\plot';
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

print(gcf, fullfile(plot_dir, 'data_analysis.png'), '-dpng', '-r300');
print(gcf, fullfile(plot_dir, 'data_analysis.jpg'), '-djpeg', '-r300');

fprintf('\n分析结果已保存到: %s\n', plot_dir);
fprintf('  - data_analysis.png\n');
fprintf('  - data_analysis.jpg\n');

% 结论
fprintf('\n=== 分析结论 ===\n');
fprintf('1. Medici数据和Matlab数据在lgN值上存在系统性差异\n');
fprintf('2. 差异随BV的增加而变化，呈现非线性关系\n');
fprintf('3. 需要采用非均匀拉伸方法来调整Medici数据，使其与Matlab数据更接近\n');
fprintf('4. 推荐使用幂函数或指数函数等非线性变换方法\n');
