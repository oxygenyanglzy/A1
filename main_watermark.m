% 1. 加载图像
original_image = imread('house.jpg'); % 替换为你的图像文件
original_image = rgb2gray(original_image); % 转换为灰度图像


N = 256; % 示例变换大小
n = 0; % 基础变换阶数

% 2. 应用离散汉克尔变换（假设有 dht 函数）
h_function = @(r) original_image; % 用一个示例函数句柄
R = max(size(original_image)); % 示例半径
H = dht(h_function, R); % 对图像进行 DHT

% 3. 逆变换（假设有 idht 函数）
recovered_image = idht(H); % 恢复图像

% 4. 显示图像
subplot(1, 2, 1);
imshow(original_image);
title('Original Image');

subplot(1, 2, 2);
imshow(recovered_image, []);
title('Recovered Image');

% 计算恢复误差
error = norm(double(original_image) - double(recovered_image), 'fro');
disp(['Recovery Error: ', num2str(error)]);

