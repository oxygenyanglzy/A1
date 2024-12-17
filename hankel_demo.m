
    % 示例信号
    signal =[159,120,57,79;
        150,59,43,236;
        53,216,58,110;
        77,49,111,47];

    % 执行 Hankel 变换
    H = hankel_transform(signal);
    disp('Hankel Matrix:');
    disp(H);

    % 执行逆 Hankel 变换
    recovered_signal = hankel_inverse_transform(H);
    disp('Recovered Signal:');
    disp(recovered_signal);
  
function H = hankel_transform(signal)
    % 计算 Hankel 矩阵
    n = length(signal);
    half = n / 2;
    % 构造 Hankel 矩阵
    H = zeros(half, half);
    for i = 1:half
        H(i, :) = signal(i:i+half-1);
    end
end

function signal = hankel_inverse_transform(H)
    % 计算逆 Hankel 变换
    [rows, cols] = size(H);
    if rows ~= cols
        error('Input matrix must be square.');
    end
    % 从 Hankel 矩阵恢复信号
    signal = zeros(1, rows + cols - 1);
    signal(1:rows) = H(:, 1)';
    signal(rows+1:end) = H(rows, 2:end);
end
