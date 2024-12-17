function psnrvalue = PSNR(original,unzip)
%% 计算原始图像的信号功率
A = imread(original);
A = double(A);
B = imread(unzip);
B = double(B);

%% 计算MSE

[m,n] = size(A);
[m1,n1] = size(B);
if m~=m1||n~=n1
    error('图像大小不一致');
end

msevalue = 0;
for i = 1:m
    for j = 1:n
        msevalue = msevalue+(A(i,j)-B(i,j))^2;
    end
end
msevalue = msevalue/(m*n);
if msevalue == 0
    error('图像完全相同');
end

%% 计算峰值信噪比
psnrvalue = 255^2/msevalue;
psnrvalue = 10*log10(psnrvalue);
