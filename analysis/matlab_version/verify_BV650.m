% 验证BV=650V时计算的W和N

% 常量定义
global a_n b_n a_p b_p q eps0 eps_r eps_s;
a_n = 7.03e5;        % cm^-1
b_n = 1.231e6;       % V/cm
a_p = 1.582e6;       % cm^-1
b_p = 2.036e6;       % V/cm
q = 1.602e-19;       % C
eps0 = 8.854e-14;    % F/cm
eps_r = 11.7;
eps_s = eps0 * eps_r;

% 计算电子电离系数
function an = alpha_n(E)
global a_n b_n;
    E = max(E, 1e-12);
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

% 计算BV=650V对应的W和N
BV = 650;
grid_points = 20000;

% 使用Fulop近似作为初始猜测
W_fulop_um = 0.0257 * (BV^(7.0/6.0));
W_guess_cm = W_fulop_um * 1e-4;
Wmin = max(1e-6, W_guess_cm * 0.2);
Wmax = max(Wmin*10, W_guess_cm * 4.0);

try
    W_cm = solve_W_for_BV(BV, Wmin, Wmax, grid_points);
catch
    W_cm = solve_W_for_BV(BV, Wmin*0.05, Wmax*20, grid_points);
end

% 计算N
N_val = 2.0 * eps_s * BV / (q * W_cm^2);

fprintf('========================================\n');
fprintf('BV = %d V 的计算结果:\n', BV);
fprintf('========================================\n');
fprintf('W = %.4f um (%.6f cm)\n', W_cm*1e4, W_cm);
fprintf('N = %.4e cm^-3\n', N_val);
fprintf('\n');
fprintf('用户给出的参考值:\n');
fprintf('W = 55.2327 um\n');
fprintf('N = 2.7557e14 cm^-3\n');
fprintf('\n');
fprintf('差异:\n');
fprintf('W 差异 = %.4f um (%.2f%%)\n', (W_cm*1e4 - 55.2327), abs(W_cm*1e4 - 55.2327)/55.2327*100);
fprintf('N 差异 = %.4e cm^-3 (%.2f%%)\n', abs(N_val - 2.7557e14), abs(N_val - 2.7557e14)/2.7557e14*100);
fprintf('\n');

% 验证电离积分是否为1
I_check = ionization_integral_for_W(W_cm, BV, grid_points);
fprintf('验证: 电离积分 I_n = %.6f (应该接近1)\n', I_check);
