function  M=Hundun( image,crypt)
%crypt=1，表示加密，cyrpt=2,
% A_encry=gongg();
%   此处显示详细说明
initial_value = [0.001,0.002,9.69];
[~,~,A] = Lorenz(initial_value);
% 由于起始时刻的数据都比较小，因此剔除起始时刻的数据
A = A(100:end);
% 归一化混沌序列
A_normal = (A-min(A))./(max(A)-min(A));
% 控制范围为[0,255]
A_uint8 = uint8(A_normal*255);
A_encry = reshape(A_uint8(1:32*32),32,32);
if  crypt==1
% 执行图像加密
% subplot(1,3,1);
% imshow(image);
 for i = 1:32
     for j = 1:32
         image_encry(i,j) = bitxor(image(i,j),A_encry(i,j));
     end
 end
        M=image_encry;

% subplot(1,3,2);
% imshow(image_encry);
% title('加密图像');

%执行图像解密
else
    image_encry=image;
for i = 1:32
    for j = 1:32
        image_decry(i,j) = bitxor(image_encry(i,j),A_encry(i,j));
    end
end
  M=image_decry;     
% subplot(1,3,3);
% imshow(image_decry);
% title('解密图像');
  end
end

