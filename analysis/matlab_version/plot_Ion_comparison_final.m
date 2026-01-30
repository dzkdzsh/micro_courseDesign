% 绘制电离积分Ion(V)对比图 - 黑白风格最终版
% 对比：Medici仿真结果 vs MATLAB计算结果 (W=55.2327μm, N=2.7557e14, BV=650V)
% 风格与R_onsp_400_1200.png保持一致：黑白、无网格、简洁图例

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

%% 器件参数（BV=650V解得的参数）
W_cm = 55.2327e-4;   % 55.2327 μm
N_dop = 2.7557e14;   % cm^-3
BV = 650;            % V
grid_points = 20000;

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

%% 计算给定外加电压V下的电离积分
function I = ionization_integral(V, W_cm, BV, grid_points)
if V <= 0
    I = 0;
    return;
end

x = linspace(0, W_cm, grid_points)';
E = 2.0 * V * (W_cm - x) / (W_cm^2);
an = alpha_n(E); an = an(:);
ap = alpha_p(E); ap = ap(:);
diff = an - ap; diff = diff(:);
S = cumtrapz(x, diff); S = S(:);
if length(S) < length(x), S = [0; S]; end
S = S(1:length(an));
integrand = an .* exp(-S);
I = trapz(x, integrand);
end

%% 读取Medici仿真数据
excel_file = 'd:\micro_courseDesign\仿真\设定650_步长0.2_击穿电压646.1\设定650_步长0.2_击穿电压646.1.xlsx';
medici_data = readtable(excel_file);
V_medici = medici_data.('VoltageV_K__V_');
Ion_medici = medici_data.('ExtractedExpression_testn_');

%% 计算MATLAB电离积分
V_calc = linspace(0, max(V_medici), 200)';
Ion_matlab = zeros(size(V_calc));

fprintf('计算MATLAB电离积分...\n');
for i = 1:length(V_calc)
    V = V_calc(i);
    Ion_matlab(i) = ionization_integral(V, W_cm, BV, grid_points);
end

fprintf('在V=650V时, Ion=%.4f (应该接近1)\n', ionization_integral(650, W_cm, BV, grid_points));

%% 创建图形 - 黑白风格
figure('Position', [100, 100, 700, 550], 'Color', 'white');

% 设置坐标轴样式（与R_onsp_400_1200.png一致）
ax = gca;
ax.Box = 'on';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.LineWidth = 0.8;
ax.FontSize = 12;
ax.FontName = 'Times New Roman';

% Medici仿真数据 - 倒三角散点（黑色）
scatter(V_medici, Ion_medici, 50, 'k', 'v', 'filled', 'DisplayName', 'I_{on}(MEDICI)');
hold on;

% MATLAB计算结果 - 黑色连续曲线
plot(V_calc, Ion_matlab, 'k-', 'LineWidth', 1.2, 'DisplayName', 'I_{on}(MATLAB)');

hold off;

% 设置坐标轴标签
xlabel('Voltage V (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Ionization Integral I_{on}', 'FontSize', 14, 'FontName', 'Times New Roman');

% 设置坐标轴范围（调整比例尺，为右上角图例留出空间）
xlim([0, max(V_medici)*1.15]);  % X轴延长15%，给图例留出空间
ylim([0, max(max(Ion_medici), max(Ion_matlab))*1.25]);  % Y轴延长25%，给图例留出足够空间

% 设置图例（右上角，与R_onsp_400_1200.png一致）
legend('Location', 'northeast', 'FontSize', 11, 'FontName', 'Times New Roman');

%% 保存图片
plot_dir = 'd:\micro_courseDesign\plot';
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

% 保存为PNG (高分辨率)
print(gcf, fullfile(plot_dir, 'Ion_V_comparison_final.png'), '-dpng', '-r300');
% 保存为JPG
print(gcf, fullfile(plot_dir, 'Ion_V_comparison_final.jpg'), '-djpeg', '-r300');

fprintf('\n图片已保存到: %s\n', plot_dir);
fprintf('  - Ion_V_comparison_final.png\n');
fprintf('  - Ion_V_comparison_final.jpg\n');

%% 显示关键数据点
fprintf('\n关键数据点对比:\n');
fprintf('%-12s %-18s %-18s\n', 'V (V)', 'Medici Ion', 'MATLAB Ion');
fprintf('%s\n', repmat('-', 1, 50));

key_V = [100, 300, 500, 600, 640, 646.1, 650, 700];
for V_target = key_V
    [~, idx_med] = min(abs(V_medici - V_target));
    [~, idx_calc] = min(abs(V_calc - V_target));
    V_val = V_medici(idx_med);
    Ion_med = Ion_medici(idx_med);
    Ion_mat = Ion_matlab(idx_calc);
    fprintf('%-12.1f %-18.6f %-18.6f\n', V_val, Ion_med, Ion_mat);
end
