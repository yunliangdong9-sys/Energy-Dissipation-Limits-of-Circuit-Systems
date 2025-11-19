clc
clear
close

% 参数设置
% n = 10000;            % 模拟的步骤数
e = 0.001;             % 错误概率
% alpha = [10, 1, 1, 10]; % 狄利克雷分布的形状参数 (较大权重分配给 x = y)
idx = 1;
T=1;
k = 1;
n = 1;



% 初始化输入和输出
% p = linspace(0, 1, 101); % 输入分布作为自变量
% p = 0.6;
p_channel = [1 - e, e; e, 1 - e]; % 信道矩阵，e是错误率
% 第一阶段的先验分布
qA_x = [0.5, 0.5]; 
% X = zeros(1, n);      % 输入
% Y = zeros(1, n);      % 输出
% Y_prev = 0;           % 上一次的输出（初始为0）
I_XY = zeros(101, 4999);

for i = 1 : 4999
    P(i) = i * 0.0001;
end
dirichlet_probs = zeros(4999, 4);
% 狄利克雷分布的概率计算，作为第二阶段的先验分布
% dirichlet_params = alpha;
% probs = gamrnd(alpha, 1, 1, numel(alpha));
% dirichlet_probs = probs / sum(probs);
for i = 1 : 4999
    dirichlet_probs(i, 1) = 0.5 - P(i);
    dirichlet_probs(i, 2) = P(i);
    dirichlet_probs(i, 3) = P(i);
    dirichlet_probs(i, 4) = 0.5 - P(i);
end
% dirichlet_matrix = reshape(dirichlet_probs, [2, 2]); % 2*2矩阵
% qB_x = sum(dirichlet_matrix, 2); % 按列求和（对Y求和），X先验概率的边缘分布
qB_x = [0.5; 0.5];

for l = 0 : 0.01 : 1
    p(idx) = l;

    % 输入分布
    p_in = [p(idx), 1 - p(idx)]; 
    % 输出分布
    p_y = p_in * p_channel; 
    % 第二阶段的联合分布
    p0_xy = [(1-e)*p_in(1), e*p_in(1); e*p_in(2), (1-e)*p_in(2)];
    %第一阶段的联合分布
    p1_xy = [p_in(1) * p_y(1), p_in(1) * p_y(2), p_in(2) * p_y(1), p_in(2) * p_y(2)]; 
    
    for m = 1 : 4999
        % 计算互信息 I(X; Y)
        for i = 1:2
            for j = 1:2
                if p0_xy(i, j) > 0
                    I_XY(idx, m) = I_XY(idx, m) + p0_xy(i, j) * log(p0_xy(i, j) / (p_in(i) * p_y(j)));
                end
            end
        end
        dirichlet = [dirichlet_probs(m,1), dirichlet_probs(m,2), dirichlet_probs(m,3), dirichlet_probs(m,4)];
        EP(idx, m) = I_XY(idx, m) + relative_entropy(p_in, qA_x) + ...
            relative_entropy(p1_xy, dirichlet) - relative_entropy(p_in, qB_x.');
        
        Q(idx, m) = k*T * (EP(idx, m));
        if (m > 1 && Q(idx, m) < Q(idx, m - 1))
            m_min(idx) = m;
        end
        eta(idx, m) = I_XY(idx, m) / Q(idx, m);
    end
    % 失配熵产
    
    idx = idx + 1;
end

for i = 1 : 1 : 50
    if (Q(i, m_min(i)) >= Q(51, m_min(i)))
        n = m_min(i);
        break;
    end
end

max_mean_idx = 1;
mean_eta = zeros(1, 4999);
for i = 1 : 4999
    for j = 2 : 100
        mean_eta(i) = mean_eta(i) + eta(j, i);
    end
    if i > 1 && mean_eta(i) > mean_eta(i - 1)
        max_mean_idx = i;
    end
end

mean_eta = mean_eta / 99;

figure;
surf(p, P, Q', 'FaceAlpha', 0.8);
xlabel('Input probability');
ylabel('Error probability of the prior distribution');
zlabel('Energy dissipation (kT)');
shading flat;
colorbar;
for i = 1 : 101
    hold on;
    plot3([p(i)], [P(m_min(i))], [Q(i, m_min(i))], '.','MarkerSize',8, 'Color', [230/255, 135/255, 133/255], 'LineWidth',3);
end
for i = 1 : 4999
    hold on;
    plot3([p(51)], [P(i)], [Q(51, i)], '.','MarkerSize',8, 'Color', "#B6B6B6", 'LineWidth',3);
end
for i = 1 : 101
    hold on;
    plot3([p(i)], [P(n)], [Q(i, n)], '.','MarkerSize',8, 'Color', [248/255, 203/255, 201/255], 'LineWidth',3);
end

figure;
plot(p, Q(:, 1), 'LineWidth', 2, 'Color', [236/255, 116/255, 117/255]);
hold on;
xlabel('Input probability');
ylabel('Energy dissipation (kT)');
grid on;
hold on;
plot([p(51)], [Q(51, 1)],'*','MarkerSize', 8, 'Color', "#D95319", 'LineWidth',1.5);

figure;
plot(p, Q(:, 4999), 'LineWidth', 2, 'Color', [236/255, 116/255, 117/255]);
hold on;
xlabel('Input probability');
ylabel('Energy dissipation (kT)');
grid on;
hold on;
plot([p(51)], [Q(51, 4999)],'*','MarkerSize', 8, 'Color', "#D95319", 'LineWidth',1.5);

figure;
plot(p, Q(:, 1000), 'LineWidth', 2, 'Color', [236/255, 116/255, 117/255]);
hold on;
xlabel('Input probability');
ylabel('Energy dissipation (kT)');
grid on;
hold on;
plot([p(11) p(51) p(91)], [Q(11, 1000) Q(51, 1000) Q(91, 1000)],'*','MarkerSize', 8, 'Color', "#D95319", 'LineWidth',1.5);

figure;
surf(P(:), p, eta(:, :));
xlabel('Error probability of the prior distribution');
ylabel('Input probability');
zlabel('Information energy efficiency (bits/kT)');
shading flat;
colorbar;
for i = 1 : 4999
    hold on;
    plot3([P(i)], [p(51)], [eta(51, i)], '.','MarkerSize', 8, 'Color', "#B6B6B6", 'LineWidth',3);
end
for i = 1 : 101
    hold on;
    plot3([P(2500)], [p(i)], [eta(i, 2500)], '.', 'MarkerSize', 8, 'Color', "#B6B6B6", 'LineWidth',3);
end

figure;
xlabel('Error probability of the prior distribution');
ylabel('Input probability');
imagesc(eta);
colorbar;

figure;
plot(P, mean_eta, 'LineWidth', 2);
xlabel('Error probability of the prior distribution');
ylabel('Average information energy efficiency (bits/kT)');
grid on;
hold on;
plot(P(max_mean_idx), mean_eta(max_mean_idx), '*','MarkerSize',10, 'Color', "#D95319", 'LineWidth',1.5);