function extracted_watermark = extract_watermark_hankel(original_image, watermarked_image, alpha)
    % 转换图像为double类型
    original_image = im2double(original_image);
    watermarked_image = im2double(watermarked_image);
    
    % 获取图像尺寸
    [M, N] = size(original_image);
    
    % 计算Hankel变换（优化版）
    original_transformed = hankel_transform(original_image);
    watermarked_transformed = hankel_transform(watermarked_image);
    
    % 提取水印
    watermark_transformed = (watermarked_transformed - alpha * original_transformed) / (1 - alpha);
    
    % 逆Hankel变换
    extracted_watermark = inverse_hankel_transform(watermark_transformed);
    
    % 显示结果
    imshow(mat2gray(extracted_watermark));
    title('Extracted Watermark');
end

function H_transformed = hankel_transform(image)
    % 示例的优化Hankel变换实现
    % 这需要替换成实际的高效实现
    H_transformed = fft2(image); % 示例代码，请用实际的Hankel变换函数
end

function image = inverse_hankel_transform(H_transformed)
    % 示例的优化逆Hankel变换实现
    % 这需要替换成实际的高效实现
    image = ifft2(H_transformed); % 示例代码，请用实际的逆Hankel变换函数
end

function H = hankel_transform(I)
    % 计算Hankel变换
    [M, N] = size(I);
    H = zeros(M, N);
    for u = 1:M
        for v = 1:N
            r = sqrt((u - M/2)^2 + (v - N/2)^2);
            theta = atan2(v - N/2, u - M/2);
            H(u, v) = sum(sum(I .* exp(-1i * r * theta)));
        end
    end
end

function I = inverse_hankel_transform(H)
    % 计算逆Hankel变换
    [M, N] = size(H);
    I = zeros(M, N);
    for u = 1:M
        for v = 1:N
            r = sqrt((u - M/2)^2 + (v - N/2)^2);
            theta = atan2(v - N/2, u - M/2);
            I(u, v) = real(sum(sum(H .* exp(1i * r * theta))));
        end
    end
end