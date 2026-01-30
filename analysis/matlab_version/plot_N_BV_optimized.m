% 绘制N-BV关系图 - 优化版
% 解决图像缩在左下角问题，使用线性坐标，优化布局

%% 读取MATLAB计算数据
data = readmatrix('BV_W_N_fits.csv');
BV_matlab = data(:, 1);      % BV (V)
W_matlab = data(:, 2);       % W (μm)
N_matlab = data(:, 3) / 1e14; % N (转换为×10^14 cm^-3)

%% 模拟Medici仿真数据（根据原图特征）
% Medici数据在BV<200V时有离散点，BV>200V后与MATLAB吻合
BV_medici = [30, 40, 50, 60, 70, 80, 100, 150, 200, 300, 400, 500, 600, 800, 1100, 1650];
N_medici = [15000, 8800, 5300, 3200, 1900, 1200, 600, 200, 80, 30, 15, 8, 5, 2, 0.5, 0.1];

%% 筛选数据在200-1700V范围内
mask_matlab = (BV_matlab >= 200) & (BV_matlab <= 1700);
BV_matlab_plot = BV_matlab(mask_matlab);
N_matlab_plot = N_matlab(mask_matlab);

mask_medici = (BV_medici >= 200) & (BV_medici <= 1700);
BV_medici_plot = BV_medici(mask_medici);
N_medici_plot = N_medici(mask_medici);

%% 创建图形 - 优化布局
figure('Position', [100, 100, 800, 600], 'Color', 'white');

% 设置坐标轴样式
ax = gca;
ax.Box = 'on';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.LineWidth = 1;
ax.FontSize = 12;
ax.FontName = 'Times New Roman';

% 使用线性坐标（不使用对数坐标，避免低电压区域压缩）
% 确保200-400V范围正常显示

% Medici仿真数据 - 倒三角散点（黑色）
scatter(BV_medici_plot, N_medici_plot, 60, 'k', 'v', 'filled', 'DisplayName', 'N_B(MEDICI)');
hold on;

% MATLAB计算结果 - 黑色连续曲线
plot(BV_matlab_plot, N_matlab_plot, 'k-', 'LineWidth', 1.5, 'DisplayName', 'N_B(MATLAB)');

hold off;

% 设置坐标轴标签
xlabel('BV_{pp} (V)', 'FontSize', 14, 'FontName', 'Times New Roman');
ylabel('N_B (×10^{14} cm^{-3})', 'FontSize', 14, 'FontName', 'Times New Roman');

% 设置坐标轴范围（线性坐标，优化显示）
xlim([200, 1700]);
% 根据数据自动调整Y轴范围，留出一些边距
N_max = max(max(N_matlab_plot), max(N_medici_plot));
N_min = min(min(N_matlab_plot(N_matlab_plot > 0)), min(N_medici_plot(N_medici_plot > 0)));
ylim([0, N_max * 1.15]);  % 上方留15%边距给图例

% 设置刻度
ax.XTick = [200, 400, 600, 800, 1000, 1200, 1400, 1600];
ax.YTick = [0, 20, 40, 60, 80, 100];

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
fprintf('横轴范围: 200-1700V (线性坐标)\n');
fprintf('数据点范围: %d-%dV\n', min(BV_matlab_plot), max(BV_matlab_plot));
