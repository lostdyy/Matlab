function J = calculate_tanh_jacobian(x)
    % 输入：
    % x: 输入向量
    % 输出：
    % J: 雅克比矩阵，即函数 y = tanh(x) 的雅克比矩阵

    % 计算 tanh 函数的输出
    y = tanh(x);

    % 计算雅克比矩阵
    n = numel(x); % 输入向量的维度
    J = zeros(n, n);

    for i = 1:n
        for j = 1:n
            if i == j
                J(i, j) = 1 - y(i)^2; % tanh 导数公式
            else
                J(i, j) = 0;
            end
        end
    end
end