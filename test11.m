% 清理工作区
clear;
clc;

% 1. 生成测试图像（假设为二维矩阵）
original_image = double(randi([0, 255], 8, 8));
h = original_image(:); % 图像数据列向量

% 2. 生成水印
watermark = randi([0, 255], 4, 4); % 生成 4x4 的随机水印矩阵

h = original_image(:); % 原始图像的列向量
h = h(:);  % 将 h 转换为列向量

% 假设我们有一个最大半径值
r_max = 1.0; % 或其他合适的最大值
num_points = length(h); % 确保点数和 h 的长度一致
r = linspace(0, r_max, num_points);
disp(diff(r));



K = max(r); % 频率范围


% 3. 对原始图像进行 Hankel 变换
[H_image, k] = fht(h, K, r, 0); % 使用 Hankel 变换将图像转换到变换域
h_watermark = watermark(:);
% 4. 嵌入水印
% 将水印嵌入到变换域
% 在 Hankel 变换域中，将水印嵌入到变换系数的左上角
H_image(1:4, 1:4) = H_image(1:4, 1:4) + double(watermark);

% 5. 逆 Hankel 变换，得到嵌入水印后的图像
watermarked_image = ifht(H_image); % 使用 Hankel 逆变换将图像转换回空间域

% 6. 对嵌入水印后的图像进行 Hankel 变换
H_watermarked = fht(double(watermarked_image));

% 7. 提取水印
% 提取水印
extracted_watermark = H_watermarked(1:4, 1:4) - H_image(1:4, 1:4);

% 将提取的水印转换为整数类型并显示
extracted_watermark = round(real(extracted_watermark));

% 8. 验证可逆性
if isequal(watermark, extracted_watermark)
    fprintf('水印提取成功，可逆性验证通过。\n');
else
    fprintf('水印提取失败，可逆性验证未通过。\n');
end

% 打印结果
disp('原始图像：');
disp(original_image);
disp('嵌入水印后的图像：');
disp(watermarked_image);
disp('提取的水印：');
disp(extracted_watermark);
function [H, k] = fht(h, K, r, n)
    % 确保 r 是均匀间距的向量
    if any(abs(diff(r) - diff(r(1:2))) > 1e-10)
        error('r 必须是均匀间距的向量');
    end

    % 频率向量 k
    k = linspace(0, K, numel(r));

    % 初始化 Hankel 变换结果
    H = zeros(size(k));

    % 确保 h 是列向量
    h = h(:);

    % 检查 h 和 r 的大小
    if length(h) ~= length(r)
        error('h 必须与 r 具有相同的长度');
    end

    % 计算每个频率点的 Hankel 变换
    for i = 1:numel(k)
        bessel_vals = besselj(n, k(i) * r);
        if length(bessel_vals) ~= length(h)
            error('besselj 计算结果与 h 的长度不匹配');
        end
        H(i) = trapz(r, h .* bessel_vals);
    end
end


function h = ifht(H, k, r, n)
    N = numel(k);
    if nargin < 4 || isempty(n)
        n = 0;
    end
    if numel(n) > 1
        I = n;
    else
        a = log(k(end)/k(1))/(N-1);
        I = [k k(end)*exp(a*(1:N))];
        I = ifft(a*r(1)*I.*besselj(n,r(1)*I));
    end
    h = fft(fft(H.*k, 2*N).*I);
    if isreal(H)
        h = real(h);
    end
    h = h(1:N)./(2*pi*r);
end

