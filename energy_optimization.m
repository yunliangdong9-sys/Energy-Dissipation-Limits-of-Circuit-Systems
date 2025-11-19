clc
clear
close

% 参数设置
e = 0.001;             % 错误概率
N = 300;                 % 比特数
t_proc = 9.99 : -0.01 : 5;         % 信息处理时间
t_tran = 0.01 : 0.01 : 5;    % 信息传输时间
c = 1;                  % 时间常数
f = 1;                  % 单位时间使用信道次数
n_max = 1000000;             % 最大信道数
Rr = N ./ t_tran;       % 信息传输速率
Rm = f * calculate_mutual_information(1/2, e); % 最大传输速率

% 初始化输入和输出
p_channel = [1 - e, e; e, 1 - e]; % 信道矩阵，e是错误率
% 第一阶段的先验分布
qA_x = [0.5, 0.5]; 
dirichlet_probs = [0.499, 0.499, 0.001, 0.001];
dirichlet_matrix = reshape(dirichlet_probs, [2, 2]); % 2*2矩阵
qB_x = sum(dirichlet_matrix, 2); % 按列求和（对Y求和），X先验概率的边缘分布

idx = 1;        % 能耗最大优化百分比索引

for i = 1 : 1 : 500
    if Rr(i) > 0
        n(i) = ceil(Rr(i)/Rm);
    end
    Q = zeros(1, n_max);
    for j = n(i) : 1 : n_max
        R(j) = Rr(i) / j;
        p(j) = find_input_probability(R(j), e);
        p_tran(i) = p(j);
        % 输入分布
        p_in = [p(j), 1 - p(j)]; 
        % 输出分布
        p_y = p_in * p_channel; 
        % 第二阶段的联合分布
        p0_xy = [(1-e)*p_in(1), e*p_in(1); e*p_in(2), (1-e)*p_in(2)];
        %第一阶段的联合分布
        p1_xy = [p_in(1) * p_y(1), p_in(1) * p_y(2), p_in(2) * p_y(1), p_in(2) * p_y(2)]; 

        I_XY = 0;
        for k = 1 : 2
            for l = 1 : 2
                if p0_xy(k, l) > 0
                    I_XY = I_XY + p0_xy(k, l) * log(p0_xy(k, l) / (p_in(k) * p_y(l)));
                end
            end
        end

        dirichlet = dirichlet_probs;
        Q(j) = j * f * t_tran(i) * (I_XY + relative_entropy(p_in, qA_x) + ...
            relative_entropy(p1_xy, dirichlet) - relative_entropy(p_in, qB_x.'));
        if j == n(i)
            Q_tran_normal(i) = Q(j);
            p_tran_normal(i) = p(j)
        end
        Q_tran(i) = Q(j);
        if j > n(i) && Q(j) > Q(j - 1) && Q(j - 1) > 0
            n(i) = j - 1;
            p_tran(i) = p(j - 1);
            Q_tran(i) = Q(j - 1);
            break;
        end
    end
    % 确定信息处理最优并行数
    M1 = ceil(N * sqrt(c / (c + p_tran(i) * t_proc(i) * log(3)))); 
    M2 = floor(N * sqrt(c / (c + p_tran(i) * t_proc(i) * log(3))));
    Q_proc1 = (N + M1) * p_tran(i) * log(3) + c * (N + M1)^2 / (M1 * t_proc(i));
    Q_proc2 = (N + M2) * p_tran(i) * log(3) + c * (N + M2)^2 / (M2 * t_proc(i));
    % 信息处理能耗
    if Q_proc1 < Q_proc2
        Q_proc(i) = Q_proc1;
    else
        Q_proc(i) = Q_proc2;
    end 
    Q_proc_normal_1(i) = N * (p_tran_normal(i) * log(3) + c * N / t_proc(i));
    % 优化前的总能耗
    Q_total_normal_1(i) = Q_proc_normal_1(i) + Q_tran_normal(i);  
    % 优化后的总能耗
    Q_total(i) = Q_proc(i) + Q_tran(i);  
    %能耗优化百分比
    percent(i) = (Q_total_normal_1(i) - Q_total(i)) / Q_total_normal_1(i);
    if percent(i) > percent(idx)
        idx = i;
    end
end

max_percent = percent(idx);

figure;
plot(t_tran, Q_total, 'LineWidth', 2, 'Color', [85/255, 146/255, 230/255]);
hold on;
plot(t_tran, Q_total_normal_1, 'LineWidth', 2, 'Color', [236/255, 116/255, 117/255]);
hold on;
xlabel('Percentage of information transmission time');
ylabel('Total energy dissipation (kT)');
legend('Optimized energy consumption', 'Energy consumption before optimization');
grid on;

function q = find_input_probability(R, e)
    % R: 已知的传输速率
    % p: 信道错误概率
    % q: 输入概率 P(X=1)

    % 定义 q 的范围
    q_range = linspace(0.5, 1, 10000);

    % 计算互信息
    I_values = arrayfun(@(q) calculate_mutual_information(q, e), q_range);

    % 找到最接近已知互信息的 q
    [~, idx] = min(abs(I_values - R));
    q = q_range(idx);
end

function I_XY = calculate_mutual_information(p, e)
    % 计算输出概率
    P_Y1 = p * (1 - e) + (1 - p) * e;
    P_Y0 = p * e + (1 - p) * (1 - e);

    % 计算输出熵 H(Y)
    H_Y = -P_Y0 * log(P_Y0) - P_Y1 * log(P_Y1);

    % 计算条件熵 H(Y|X)
    H_Y_given_X = -e * log(e) - (1 - e) * log(1 - e);

    % 计算互信息
    I_XY = H_Y - H_Y_given_X;
end
