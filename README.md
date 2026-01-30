# micro_courseDesign

**推导与数值验证（针对论文New expressions for non-pu...ons based on Chynoweth law_黄海猛 Eq.(3)-(5)）**

**目的：** 比较论文中 N_pp(BV)、W_pp(BV) 与 W_pp(N_pp) 的经验式，并给出推导要点、数值流程与对比结果。

**1 断裂积分与 Chynoweth 模型**

论文断裂条件（式(1) 原式）采用离子化积分：

$$I_n=\int_0^W \alpha_n(E(x))\exp\Big[-\int_0^x(\alpha_n-\alpha_p)\,dv\Big]\,dx=1. $$

其中 Chynoweth 形式：

$$\alpha_n=a_n\exp(-b_n/E),\quad \alpha_p=a_p\exp(-b_p/E).$$

参数采用论文值： $a_n=7.03\times10^5\ \mathrm{cm^{-1}}$，$b_n=1.231\times10^6\ \mathrm{V/cm}$，$a_p=1.582\times10^6\ \mathrm{cm^{-1}}$，$b_p=2.036\times10^6\ \mathrm{V/cm}$。


**2 NPT 情形的电场关系（用于消去 N）**

NPT 中线性场：

$$E(x)=\frac{qN(W-x)}{\varepsilon_s}=\frac{2\,BV\,(W-x)}{W^2},$$

并由此得：

$$BV=\frac{q N W^2}{2\varepsilon_s},\quad N=\frac{2\varepsilon_s\,BV}{qW^2}.$$ 

于是对给定 $BV$，我们只需在 $W$ 上解方程 $I_n(W;BV)=1$，解得 $W$ 后通过上式计算 $N$，对一系列 $BV$ 做幂律拟合即得论文的经验式。


**3 数值实现要点（在 `analysis/matlab_version` 目录中的实现）**

- 代码文件： [matlab_version/compute_fits.m](matlab_version/compute_fits.m)、[matlab_version/complete_test.m](matlab_version/complete_test.m)。
- 积分方法：累积梯形积分（MATLAB 的 `cumtrapz` 与 `trapz`）。
- 求根方法：fzero 法（MATLAB `fzero`），容限采用默认值。 
- 空间网格：为保证精度采用 20000 个网格点。
- BV 采样：默认 400–1600 V，均匀采样 41 个点；为测试拟合对区间敏感性也做了 400–800、400–1000、600–1600 等子区间拟合。

运行命令（在 `analysis/matlab_version` 目录下）：

```powershell
matlab -batch "cd('d:\micro_courseDesign\analysis\matlab_version'); run('compute_fits.m'); exit;"
```


**4 拟合流程**

1. 对每个 BV，在 $W$ 区间（初始用 Fulop 近似给出搜索区间）内求解 $I_n(W;BV)=1$ 得到 $W$（单位 cm），转为 μm。
2. 计算 $N=2\varepsilon_s BV/(qW^2)$（单位 cm^-3）。
3. 对数据点做对数线性拟合：$\log N = \log C + p\log BV$，得到 $C$ 与 $p$（论文的 Eq.(3) 就是此形式）。


**5 计算结果与论文对比**

- 论文给出：
  - (3) $N_{pp}\approx 1.202\times10^{18}\,BV^{-1.292}\ \mathrm{cm^{-3}}$
  - (4) $W_{pp}\approx 0.03319\,BV^{1.146}\ \mathrm{\mu m}$
  - (5) $W_{pp}\approx 4.322\times10^{14}\,N_{pp}^{-0.8927}\ \mathrm{\mu m}$

- 在 `analysis/matlab_version` 目录中通过严格数值流程得到的拟合结果为：
  - 全区间 400–1600 V（密集点、网格20000）：

    $$N\approx 1.4316\times10^{18}\,BV^{-1.32091}\ (p=-1.320914)$$

  - 区间 600–1600 V 拟合：

    $$N\approx 1.3436\times10^{18}\,BV^{-1.31192}\ (p=-1.311919)$$

  - 对应的 W 与 N-W 的拟合結果為：

    $$W\approx 0.030056\,BV^{1.16046}\ \mathrm{\mu m},$$

    $$W\approx 2.6794\times10^{14}\,N^{-0.87851}\ \mathrm{\mu m}.$$ 

  - 不同区间的 N-BV 拟合对比：

    | 区间 (V) | 系数 C |
    | -------- | ------ ||  || --------- |
    | 400-1600 | 1.4316e+18 | -1.320914 |
    | 400-1200 | 1.5130e+18 | -1.329445 |
    | 400-1000 | 1.5591e+18 | -1.334161 |
    | 500-1600 | 1.3761e+18 | -1.315297 |
    | 600-1600 | 1.3436e+18 | -1.311919 |
    | 400-800  | 1.6273e+18 | -1.340971 ||

**6 MATLAB 代码核心逻辑**

```matlab
% 计算电子电离系数
function an = alpha_n(E)
global a_n b_n;
    E = max(E, 1e-12);  % 避免除以零
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

% 批量计算
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
```

**7 复现步骤**

1. 打开 MATLAB 软件，切换到目录 `d:\micro_courseDesign\analysis\matlab_version`
2. 运行完整计算和拟合脚本：
   ```matlab
   run('compute_fits.m');
   ```
3. 查看命令窗口输出的拟合结果
4. 检查生成的数据文件 `BV_W_N_fits.csv`

或使用 MATLAB 批处理模式（无需打开 GUI）：
```powershell
cd d:\micro_courseDesign\analysis\matlab_version
matlab -batch "run('compute_fits.m'); exit;"
```

**8 与 Python 版本的差异说明**

MATLAB 实现与 Python 版本在数值方法上略有不同：
- 积分方法：MATLAB 使用 `cumtrapz` 和 `trapz`，Python 使用 SciPy 的 `cumulative_trapezoid` 和 `trapezoid`
- 求根方法：MATLAB 使用 `fzero`，Python 使用 `scipy.optimize.brentq`
- 拟合结果：两者计算结果非常接近，差异小于 1%，主要源于数值方法的实现细节不同

MATLAB 版本的优势：
- 无需额外安装第三方库，MATLAB 自带所有需要的函数
- 代码结构清晰，易于理解和修改
- 适合在 MATLAB 环境中进行进一步的数据分析和可视化

---

**9 导通电阻 R_on,pp(BV) 的计算**

**9.1 数学模型**

在功率器件设计中，导通电阻是一个关键参数。对于非穿通（NPT）结构，比导通电阻（specific on-resistance）可以表示为：

$$R_{on,pp} = \frac{W_{pp}}{q \cdot \mu_n \cdot N_{pp}}$$

其中：
- $R_{on,pp}$：比导通电阻，单位 Ω·cm²
- $W_{pp}$：耗尽区宽度，单位 cm
- $q$：电子电荷，$1.602\times10^{-19}$ C
- $\mu_n$：电子迁移率，单位 cm²/(V·s)
- $N_{pp}$：掺杂浓度，单位 cm⁻³

**电子迁移率模型（Caughey-Thomas 模型）[14, 15]**

电子迁移率与掺杂浓度相关，采用 Caughey-Thomas 经验模型：

$$\mu_n = \mu_{n0} + \frac{\mu_{nmax}}{1 + \left(\frac{N}{N_{ref}}\right)^{\alpha}}$$

其中参数值为：
- $\mu_{n0} = 55.24$ cm²/(V·s)：低掺杂极限迁移率
- $\mu_{nmax} = 1373.99$ cm²/(V·s)：迁移率修正幅值
- $N_{ref} = 1.072\times10^{17}$ cm⁻³：参考掺杂浓度
- $\alpha = 0.73$：幂指数

**参考文献：**
- [14] Caughey, D. M., & Thomas, R. E. (1967). Carrier mobilities in silicon empirically related to doping and field. Proceedings of the IEEE, 55(12), 2192-2193.
- [15] Arora, N. D., Hauser, J. R., & Roulston, D. J. (1982). Electron and hole mobilities in silicon as a function of concentration and temperature. IEEE Transactions on Electron Devices, 29(2), 292-295.

**9.2 实现方法**

R_on,pp(BV) 的计算流程：

1. **计算 W(BV) 和 N(BV)**：使用已有的电离积分方法，在 400-1600 V 范围内计算每个 BV 对应的 W 和 N
2. **计算电子迁移率**：使用 Caughey-Thomas 模型，根据 N 计算 μ_n
3. **计算导通电阻**：代入 R_on,pp = W / (q · μ_n · N)
4. **幂律拟合**：对 R_on,pp(BV) 进行对数线性拟合，得到经验公式

**9.3 计算结果**

- 代码文件：[matlab_version/compute_Ron_pp.m](matlab_version/compute_Ron_pp.m)

运行命令：
```powershell
matlab -batch "cd('d:\micro_courseDesign\analysis\matlab_version'); run('compute_Ron_pp.m'); exit;"
```

**拟合结果（全区间 400-1600 V）：**

$$R_{on,pp} \approx 9.8848\times10^{-9} \times BV^{2.471753}\ \Omega\cdot\mathrm{cm}^2$$

其中指数 $p = 2.471753$。

**部分计算数据：**

| BV (V) | W (μm) | N (cm⁻³) | μ_n (cm²/Vs) | R_on,pp (Ω·cm²) |
|--------|--------|----------|--------------|-----------------|
| 400 | 31.284 | 5.286×10¹⁴ | 1401.37 | 0.0264 |
| 700 | 60.292 | 2.490×10¹⁴ | 1413.01 | 0.1070 |
| 1000 | 91.201 | 1.555×10¹⁴ | 1417.69 | 0.2583 |
| 1300 | 123.423 | 1.104×10¹⁴ | 1420.23 | 0.4915 |
| 1600 | 156.651 | 8.432×10¹³ | 1421.82 | 0.8156 |

**不同区间的拟合对比：**

| 区间 (V) | 系数 C | 指数 p |
|----------|--------|--------|
| 400-1600 | 9.8848×10⁻⁹ | 2.471753 |
| 400-1200 | 9.1945×10⁻⁹ | 2.482925 |
| 400-1000 | 8.8449×10⁻⁹ | 2.489010 |
| 500-1600 | 1.0403×10⁻⁸ | 2.464488 |
| 600-1600 | 1.0735×10⁻⁸ | 2.460057 |
| 400-800 | 8.3763×10⁻⁹ | 2.497676 |

**9.4 结果分析**

1. **指数特性**：R_on,pp 与 BV 的幂律关系指数约为 2.47，这意味着击穿电压每增加一倍，导通电阻约增加 2^(2.47) ≈ 5.5 倍。这反映了击穿电压与导通电阻之间的基本权衡关系。

2. **迁移率变化**：在计算范围内，电子迁移率从约 1401 cm²/(V·s) 变化到 1422 cm²/(V·s)，变化幅度较小（约 1.5%）。这是因为掺杂浓度范围（8.4×10¹³ ~ 5.3×10¹⁴ cm⁻³）处于 Caughey-Thomas 模型的高迁移率平台区。

3. **区间敏感性**：不同 BV 区间的拟合指数略有差异，范围在 2.46-2.50 之间。低 BV 区间（400-800 V）的指数略高，高 BV 区间（600-1600 V）的指数略低。

4. **物理意义**：R_on,pp ∝ BV^(2.47) 的关系比理想的 2.0-2.5 范围略高，这主要是由于 W ∝ BV^1.16 和 N ∝ BV^(-1.32) 的非理想幂律特性导致的。

**9.5 MATLAB 代码实现**

```matlab
% 计算电子迁移率（Caughey-Thomas 模型）
function mu_n = electron_mobility_ct(N)
global mu_n0 mu_nmax N_ref alpha_CT;
    N = max(N, 1e-10);
    mu_n = mu_n0 + mu_nmax ./ (1 + (N ./ N_ref).^alpha_CT);
end

% 计算导通电阻 R_on,pp
function R_on = calculate_Ron_pp(W_cm, N)
global q;
    mu_n = electron_mobility_ct(N);
    R_on = W_cm ./ (q .* mu_n .* N);
end
```

**9.6 验证与讨论**

1. **单位一致性验证**：
   - W：cm
   - q：C
   - μ_n：cm²/(V·s)
   - N：cm⁻³
   - R_on = cm / (C · cm²/(V·s) · cm⁻³) = cm / (C/(V·s·cm)) = Ω·cm²
   单位验证正确。

2. **数值合理性**：
   - 计算得到的 R_on,pp 值在 0.026-0.82 Ω·cm² 范围内，与文献报道的硅功率器件典型值一致。
   - 迁移率值在 1400-1422 cm²/(V·s) 范围内，符合低掺杂 n 型硅的电子迁移率特性。

3. **模型局限性**：
   - Caughey-Thomas 模型假设室温（300 K）条件，未考虑温度效应。
   - 模型适用于磷掺杂 n 型硅，对于其他掺杂类型（如砷、锑）可能需要调整参数。
   - 高电场下的速度饱和效应未在模型中考虑。

