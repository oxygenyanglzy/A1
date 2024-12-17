function [h_watermarked] = embed_watermark(h, watermark, K, r, n, alpha)
    % 计算原始函数的Hankel变换
    [H, k] = fht(h, K, r, n);

    % 确保 watermark 的长度与 H 一致
    if length(watermark) ~= length(H)
        error('水印的长度与原始信号的Hankel变换长度不匹配');
    end
    
    % 嵌入水印
    H_watermarked = H + alpha * watermark;
    
    % 计算逆Hankel变换以获得水印后的信号
    h_watermarked = ifht(H_watermarked, K, r, n);
end
