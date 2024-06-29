function J = calculate_tanh_jacobian(x)
    % ���룺
    % x: ��������
    % �����
    % J: �ſ˱Ⱦ��󣬼����� y = tanh(x) ���ſ˱Ⱦ���

    % ���� tanh ���������
    y = tanh(x);

    % �����ſ˱Ⱦ���
    n = numel(x); % ����������ά��
    J = zeros(n, n);

    for i = 1:n
        for j = 1:n
            if i == j
                J(i, j) = 1 - y(i)^2; % tanh ������ʽ
            else
                J(i, j) = 0;
            end
        end
    end
end