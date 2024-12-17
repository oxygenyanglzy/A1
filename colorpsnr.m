function psnrvalue = colorpsnr(original,unzip)
%% 计算原始图像的信号功率
A = double(original);
B = double(unzip);

%% 计算MSE

[m,n] = size(A);
[m1,n1] = size(B);
if m~=m1||n~=n1
    error('图像大小不一致');
end
sum=0;
for level=1:3
    msevalue = 0;
    for i = 1:m
        for j = 1:m
            msevalue = msevalue+(A(i,j,level)-B(i,j,level))^2;
        end
    end
    msevalue = msevalue/(m*m);
    if msevalue == 0
        error('图像完全相同');
    end

    %% 计算峰值信噪比
    value = 255^2/msevalue;
    value = 10*log10(value);
sum=sum+value;
end
psnrvalue=sum/3;
