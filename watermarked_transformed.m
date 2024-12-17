function [H_image, k] = watermarked_transformed(original_image, K)
    % 将图像数据展开为列向量
    h = original_image(:); % 图像数据列向量
    
    % 生成均匀间距的半径向量
    r_max = 1.0; % 你可以根据需要设置
    num_points = length(h); % 确保 r 的长度与 h 匹配
    r = linspace(0, r_max, num_points); % 生成均匀间距的 r

    % 确保 r 是均匀间距的
    if any(diff(r) ~= diff(r(1:2)))
        error('r 必须是均匀间距的向量');
    end

    % 调用 fht 函数
    [H_image, k] = fht(h, K, r, 0); % 使用 Hankel 变换将图像转换到变换域
end