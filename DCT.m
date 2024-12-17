% 读取原始图像和水印
original_image = imread('house.jpg'); % 原始图像
watermark_image = imread('ldu3232.jpg'); % 水印图像

% 转换为灰度图像
original_image = rgb2gray(original_image);
watermark_image = rgb2gray(watermark_image);

% 确保水印大小不超过原始图像
[rows_img, cols_img] = size(original_image);
[rows_wm, cols_wm] = size(watermark_image);
if rows_wm > rows_img || cols_wm > cols_img
    error('Watermark size exceeds original image size.');
end

% DCT 变换
dct_image = dct2(original_image);

% 嵌入水印（使用DCT系数的低频部分）
dct_watermarked = dct_image;
for i = 1:rows_wm
    for j = 1:cols_wm
        dct_watermarked(i, j) = dct_image(i, j) + 0.1 * double(watermark_image(i, j)); % 水印强度
    end
end

% 反 DCT 变换
watermarked_image = idct2(dct_watermarked);

% 显示结果
figure;
subplot(1, 3, 1), imshow(original_image), title('Original Image');
subplot(1, 3, 2), imshow(uint8(watermarked_image)), title('Watermarked Image');

% 保存水印图像
imwrite(uint8(watermarked_image), 'watermarked_image.png');
% 提取水印
dct_watermarked_extracted = dct2(uint8(watermarked_image));
extracted_watermark = zeros(rows_wm, cols_wm);

for i = 1:rows_wm
    for j = 1:cols_wm
        extracted_watermark(i, j) = (dct_watermarked_extracted(i, j) - dct_image(i, j)) / 0.1; % 水印强度
    end
end

% 显示提取的水印
figure;
imshow(uint8(extracted_watermark)), title('Extracted Watermark');
