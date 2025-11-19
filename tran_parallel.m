clc
clear
close

% 参数设置
f = 1;                  % 单位时间信道使用次数
e = 0.001;             % 错误概率
t = 1;                  % 传输时间

% 初始化输入和输出
p_channel = [1 - e, e; e, 1 - e]; % 信道矩阵，e是错误率
% 第一阶段的先验分布
qA_x = [0.5, 0.5]; 
dirichlet_probs = [0.499, 0.499, 0.001, 0.001];
dirichlet_matrix = reshape(dirichlet_probs, [2, 2]); % 2*2矩阵
qB_x = sum(dirichlet_matrix, 2); % 按列求和（对Y求和），X先验概率的边缘分布

Rr = 0 : 0.01 : 25;               % 目标传输速率
n = ones(1, 2501);               % 对应信道数
Rm = f * calculate_mutual_information(1/2, e);   % 最大传输速率

n_max = 70;

p_max = 0;
p_min = 0;

for i = 1 : 2501
    if Rr(i) > 0
        n(i) = ceil(Rr(i)/Rm);
    end
    for j = n(i) : 1 : n_max
        R(j) = Rr(i) / j;
        p(j) = find_input_probability(R(j), e);
        p_opt(i) = p(j);
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
        Q(j) = j * f * t * (I_XY + relative_entropy(p_in, qA_x) + ...
            relative_entropy(p1_xy, dirichlet) - relative_entropy(p_in, qB_x.'));
        Qt(i) = Q(j);
        if j > n(i) && Q(j) > Q(j - 1)
            n(i) = j - 1;
            p_opt(i) = p(j - 1);
            Qt(i) = Q(j - 1);
            if (i > 2 && p_opt(i) < p_opt(i - 1) && p_opt(i - 1) > p_opt(i - 2))
                p_max = p_opt(i - 1);
            end
            if (i > 2 && p_opt(i) > p_opt(i - 1) && p_opt(i - 1) < p_opt(i - 2))
                p_min = p_opt(i - 1);
            end
            break;
        end
    end
end

p_lim = (p_max + p_min) / 2;

figure;
plot(Rr(1:701), n(1:701), 'LineWidth', 2, 'Color', [69/255, 105/255, 144/255]);
xlabel('Target transmission rate (bits/s)');
ylabel('Optimal number of channels');
grid on;

figure;
plot(Rr, p_opt, 'LineWidth', 2, 'Color', [120/255, 183/255, 201/255]);
xlabel('Target transmission rate (bits/s)');
ylabel('Optimal input probability');
grid on;
hold on;
yline(p_lim, 'LineWidth', 2 ,"Color",[229/255, 139/255, 123/255]);

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