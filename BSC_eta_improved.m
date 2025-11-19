clc
clear
close

% 参数设置
e = 0.001;             % 错误概率
idx = 1;

% 初始化输入和输出
p_channel = [1 - e, e; e, 1 - e]; % 信道矩阵，e是错误率
% 第一阶段的先验分布
qA_x = [0.5, 0.5]; 
dirichlet_probs = [0.499999, 0.000001, 0.000001, 0.499999];
dirichlet_matrix = reshape(dirichlet_probs, [2, 2]); % 2*2矩阵
qB_x = sum(dirichlet_matrix, 2); % 按列求和（对Y求和），X先验概率的边缘分布
a = 1;
c = 1;
for i = 0 : 0.001 : 1
    p(idx) = i;
    % 输入分布
    p_in = [p(idx), 1 - p(idx)]; 
    % 输出分布
    p_y = p_in * p_channel; 
    % 第二阶段的联合分布
    p0_xy = [(1-e)*p_in(1), e*p_in(1); e*p_in(2), (1-e)*p_in(2)];
    %第一阶段的联合分布
    p1_xy = [p_in(1) * p_y(1), p_in(1) * p_y(2), p_in(2) * p_y(1), p_in(2) * p_y(2)]; 
    
    I_XY(idx) = 0;
    for k = 1 : 2
        for j = 1 : 2
            if p0_xy(k, j) > 0
                I_XY(idx) = I_XY(idx) + p0_xy(k, j) * log(p0_xy(k, j) / (p_in(k) * p_y(j)));
            end
        end
    end

    dirichlet = dirichlet_probs;
    EP(idx) = I_XY(idx) + relative_entropy(p_in, qA_x) + ...
            relative_entropy(p1_xy, dirichlet) - relative_entropy(p_in, qB_x.');
        
    Q(idx) = EP(idx);
    eta(idx) = I_XY(idx) / Q(idx);
    if (idx > 1 && idx < 500 && eta(idx) > eta(idx - 1))
        a = idx;
    end
    if p(idx) == 1 - p(a)
        c = idx;
    end
    
    idx = idx + 1;
end
    % 失配熵产
b = (eta(a) - eta(501)) / eta(501);


figure;
plot(p, eta, 'LineWidth', 2);
xlabel('Input probability');
ylabel('Information energy efficiency (bits/kT)');
grid on;
hold on;
plot([p(a), p(501)], [eta(a) eta(501)],'*','MarkerSize',8, 'Color', "#D95319", 'LineWidth',1.5);

