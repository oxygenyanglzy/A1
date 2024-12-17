function A_encry= gongg(  )
%GONGG 此处显示有关此函数的摘要
%   此处显示详细说明
initial_value = [0.001,0.005,0.002];
[~,~,A] = Lorenz(initial_value);
% 由于起始时刻的数据都比较小，因此剔除起始时刻的数据
A = A(100:end);
% 归一化混沌序列
A_normal = (A-min(A))./(max(A)-min(A));
% 控制范围为[0,255]
A_uint8 = uint8(A_normal*255);
A_encry = reshape(A_uint8(1:32*32),32,32);
end

