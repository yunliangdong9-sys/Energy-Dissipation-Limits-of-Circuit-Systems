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
dirichlet_probs = [0.25, 0.25, 0.25, 0.25];
dirichlet_matrix = reshape(dirichlet_probs, [2, 2]); % 2*2矩阵
qB_x = sum(dirichlet_matrix, 2); % 按列求和（对Y求和），X先验概率的边缘分布

idx_min = idx;
idx_max = idx;

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
    E_t(idx) = I_XY(idx) + relative_entropy(p_in, qA_x) + ...
            relative_entropy(p1_xy, dirichlet) - relative_entropy(p_in, qB_x.');
        
    E_p(idx) = log(3) * p(idx)*I_XY(idx);
    Q(idx) = E_t(idx) + E_p(idx);
    if Q(idx) < Q(idx_min) 
        idx_min = idx;
    end
    if Q(idx) > Q(idx_max)
        idx_max = idx;
    end
    eta(idx) = I_XY(idx) / Q(idx);
    
    idx = idx + 1;
end

percent_max = (Q(idx_max) - Q(idx_min)) / Q(idx_max);
percent = (Q(501) - Q(idx_min)) / Q(501);

figure;
plot(p, Q, 'LineWidth', 2, 'Color', [120/255, 183/255, 201/255]);
xlabel('p_{out}^{proc}(1)');
ylabel('Total energy dissipation (kT)');
grid on;
hold on;
plot([p(idx_min), p(idx_max), p(501)], [Q(idx_min), Q(idx_max), Q(501)],'.','MarkerSize',15, 'Color', [229/255, 139/255, 123/255], 'LineWidth',1.5);
yline([Q(idx_min), Q(idx_max), Q(501)], '--', 'LineWidth', 1.5 ,"Color",[229/255, 99/255, 123/255]);

for i = 1 : 9
    j = i * 100 + 1;
    Et_percent(i) = E_t(j) / Q(j);
end

figure;
h = bar([1 - Et_percent; Et_percent]', "stacked");
l = legend('information processing','information transmission');
set(gca,'XTick',1:9,'XTickLabel',{'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'});
set(h(2), 'FaceColor', [233, 196, 106]/256);
set(h(1), 'FaceColor', [41, 157, 143]/256);
set(h,'edgecolor','none');
set(l, 'edgecolor', [0.8,0.8,0.8]);
xlabel('p_{out}^{proc}(1)');
ylabel('Percentage of processing/transmission energy dissipation');
grid on;
box on;


