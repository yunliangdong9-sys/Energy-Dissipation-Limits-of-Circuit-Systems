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
dirichlet_probs = [0.499, 0.499, 0.001, 0.001];
dirichlet_matrix = reshape(dirichlet_probs, [2, 2]); % 2*2矩阵
qB_x = sum(dirichlet_matrix, 2); % 按列求和（对Y求和），X先验概率的边缘分布

a = 1;
b = 1;
d = 1;
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
    if abs(Q(idx) - 1.5) < abs(Q(a) - 1.5) 
        a = idx;
    end
    if abs(Q(idx) - 3) < abs(Q(b) - 3)
        b = idx;
    end
    if abs(I_XY(idx) - 0.6852 / 2) < abs(I_XY(d) - 0.6852 / 2) && Q(idx) < 3
        d = idx;
    end
    idx = idx + 1;
end
    % 失配熵产
c = (2 * I_XY(a) - I_XY(b)) / I_XY(b);
e = (Q(b) - 2 * Q(d)) / Q(b);
figure;
plot(Q, I_XY, 'LineWidth', 2, 'Color', [120/255, 183/255, 201/255]);
xlabel('Energy dissipation (kT)');
ylabel('Mutual information (bits)');
grid on;
hold on;
plot([Q(a), Q(a), Q(b), Q(d), Q(d) * 2], [I_XY(a), 2 * I_XY(a), I_XY(b), I_XY(d), I_XY(d)],'.','MarkerSize',15, 'Color', [229/255, 139/255, 123/255], 'LineWidth',1.5);
yline([I_XY(d)], ':', 'LineWidth', 1.5 ,"Color",[229/255, 99/255, 123/255]);
xline([Q(a)], ':', 'LineWidth', 1.5 ,"Color",[229/255, 99/255, 123/255]);

