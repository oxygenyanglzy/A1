function psnrvalue = PSNR(original,unzip)
%% ����ԭʼͼ����źŹ���
A = imread(original);
A = double(A);
B = imread(unzip);
B = double(B);

%% ����MSE

[m,n] = size(A);
[m1,n1] = size(B);
if m~=m1||n~=n1
    error('ͼ���С��һ��');
end

msevalue = 0;
for i = 1:m
    for j = 1:n
        msevalue = msevalue+(A(i,j)-B(i,j))^2;
    end
end
msevalue = msevalue/(m*n);
if msevalue == 0
    error('ͼ����ȫ��ͬ');
end

%% �����ֵ�����
psnrvalue = 255^2/msevalue;
psnrvalue = 10*log10(psnrvalue);
