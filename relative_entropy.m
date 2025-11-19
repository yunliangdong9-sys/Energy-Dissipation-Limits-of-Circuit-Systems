function D_kl = relative_entropy(P, Q)
    D_kl = 0;
    for i = 1:length(P)
        if P(i) > 0  % 避免 log(0) 或无意义计算
            D_kl = D_kl + P(i) * log(P(i) / Q(i));
        end
    end
end