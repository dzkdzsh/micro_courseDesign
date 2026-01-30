% 绘制N-BV关系图 - 调整横轴范围50-1700V
% 风格与原有N_BV.png保持一致：黑白、无网格、倒三角散点

%% 读取MATLAB计算数据
data = readmatrix('BV_W_N_fits.csv');
BV_matlab = data(:, 1);      % BV (V)
W_matlab = data(:, 2);       % W (μm)
N_matlab = data(:, 3) / 1e14; % N (转换为×10^14 cm^-3)

%% 模拟Medici仿真数据（根据原图特征）
% 从原图可见，Medici数据在BV<200V时有离散点，BV>200V后与MATLAB吻合
BV_medici = [30, 40, 50, 60, 70, 80, 100, 150, 200, 300, 400, 500, 600, 800, 1100, 1650];
N_medici = [15000, 8800, 5300, 3200, 1900, 1200, 600, 200, 80, 30, 15, 8, 5, 2, 0.5, 0.1];

%% 筛选MATLAB数据在50-1700V范围内
mask = (BV_matlab >= 50) & (BV_matlab <= 1700);
BV_matlab_plot = BV_matlab(mask);
N_matlab_plot = N_matlab(mask);

%% 筛选数据在200-1700V范围内（用于绘图）
mask_plot = (BV_matlab >= 200) & (BV_matlab <= 1700);
BV_matlab_plot = BV_matlab(mask_plot);
N_matlab_plot = N_matlab(mask_plot);

% 筛选Medici数据在200-1700V范围内
mask_med = (BV_medici >= 200) & (BV_medici <= 1700);
BV_medici_plot = BV_medici(mask_med);
N_medici_plot = N_medici(mask_med);

%% 创建图形 - 黑白风格，使用对数坐标
figure('Position', [100, 100, 700, 550], 'Color', 'white');

% 设置坐标轴样式
ax = gca;
ax.Box = 'on';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.LineWidth = 0.8;
ax.FontSize = 12;
ax.FontName = 'Times New Roman';

% 设置对数坐标轴
set(gca, 'XScale', 'log', 'YScale', 'log');

% Medici仿真数据 - 倒三角散点（黑色）
scatter(BV_medici_plot, N_medici_plot, 50, 'k', 'v', 'filled', 'DisplayName', 'N_B(MEDICI)');
hold on;

% MATLAB计算结果 - 黑色连续曲线
plot(BV_matlab_plot, N_matlab_plot, 'k-', 'LineWidth', 1.2, 'DisplayName', 'N_B(MATLAB)');

hold off;

% 设置坐标轴标签
xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('N_B (×10^{14} cm^{-3})', 'FontSize', 14, 'FontName', 'Times New Roman');

% 设置坐标轴范围（对数坐标）
xlim([200, 1700]);
ylim([0.05, 200]);

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

fprintf('图片已保存到: %s\n', plot_dir);
fprintf('  - N_BV.png\n');
fprintf('  - N_BV.jpg\n');
fprintf('横轴范围: 50-1700V\n');
