function rotateimg=RotatingAttack(img,angle)
%旋转攻击图像，并校正：2019-11-24完成
I2=imrotate(img,angle);
% M=512;N=512;
% if angle==0
%     rotateimg=img;
%     disp('没有对图像进行旋转攻击!');
% else
%     [r2,c2,~]=size(I2);
%     %%旋转角度测试% 大角度旋转测试
%     figure(10);
%     subplot(2,3,1),imshow(uint8(img));title('lena');
%     for j=1:c2
%         for i=1:r2
%             if I2(i,j,1)>0 %测试旋转图像与左边缘的交点
%                 ry=i;           %在垂直方向上优先定位
%                 rx=r2-ry;
%                 break;
%             end
%         end
%     end
%     angle=round(atan(rx/ry)*180/pi);%计算旋转角度
%     
%      I3=imrotate(I2,-angle);%%%利用所求角度进行逆向旋转
%     subplot(2,3,2),imshow(uint8(I2));title(['rotated image ','angle=',num2str(angle)]);
%     subplot(2,3,3),imshow(uint8(I3));title('Inverse rotation image');
%     %%----------------裁剪边缘-----------------
%     %%-------裁剪边缘，去掉黑边并缩放到指定大小[M,N]-------%%
% [f,c]=find(I3(:,:,3));
% minf=min(f);maxf=max(f);
% minc=min(c);maxc=max(c);
% 
% %测试图片左边非零元素的个数
% count=0;
% for j=minc:maxc
%     for i=minf:maxf
%         if I3(i,j,3)>=30 && I3(i,j,2)>=30 && I3(i,j,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(maxf-minf)/2  %如果首列元素半数以上是<30的元素，则说明该列无效，列数+1，直至有效为止
%         minc=minc+1;
%     else
%         break;
%     end
% end
% %测试图片右边非零元素的个数
% count=0;
% for j=minc:maxc
%     jj=maxc+minc-j;
%     for i=minf:maxf
%         if I3(i,jj,3)>=30 && I3(i,jj,2)>=30 && I3(i,jj,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(maxf-minf)/2 %如果尾列元素半数以上是<30的元素，则说明该列无效，列数-1，直至有效为止
%         maxc=maxc-1;
%     else
%         break;
%     end
% end
% %测试图片上边非零元素的个数
% count=0;
% for i=minf:maxf
%     for j=minc:maxc
%         if I3(i,j,3)>=30 && I3(i,j,2)>=30 && I3(i,j,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(maxc-minc)/2 %如果首行元素半数以上是<30的元素，则说明该行无效，行数+1，直至有效为止
%         minf=minf+1;
%     else
%         break;
%     end
% end
% %测试图片下边非零元素的个数
% count=0;
% for i=minf:maxf
%     ii=minf+maxf-i;
%     for j=minc:maxc
%         if I3(ii,j,3)>=30 && I3(ii,j,2)>=30 && I3(ii,j,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(maxc-minc)/2 %如果尾行元素半数以上是<30的元素，则说明该行无效，行数-1，直至有效为止
%         maxf=maxf-1;
%     else
%         break;
%     end
% end
% I4=I3(minf:maxf,minc:maxc,:);
% correctimg=imresize(I4,[M N]);%%%对裁边的图像进行缩放到规则尺寸
% subplot(2,3,4),imshow(uint8(I4)),title('裁剪边缘');
% correctimg=imresize(I4,[512 512]);%%%对裁边的图像进行缩放到规则尺寸
% psnr=colorpsnr(watermarkedimg,correctimg);%%计算校准图像与原始图像Host的相似度PSNR
% disp(['测得的旋转角度：',num2str(angle),'    PSNR=',num2str(psnr)]);    %%求的旋转角度
% subplot(2,3,5),imshow(uint8(correctimg));title(['规范尺寸 PSNR=',num2str(psnr,-4)]);
% end
rotateimg=I2;
% figure;
% imshow(uint8(rotateimg));
% imwrite(uint8(I2),['Lena-Rotation',num2str(angle),'.jpg']);
end
% % % [r3,c3,~]=size(I3);
% % %     flag=0; %定位图像最左列
% % %     for j=1:c3
% % %         for i=1:r3
% % %             if I3(i,j,1)>0  %%在列的方向上找到第一个不为0 的点的坐标[rr1,cc1]
% % %                 %         rr1=i;          %在列方向上定位最左边
% % %                 cc1=j;
% % %                 flag=1;
% % %                 break;
% % %             end
% % %         end
% % %         if flag==1
% % %             break;
% % %         end
% % %     end
% % %     
% % %     flag=0;%定位图像最上行
% % %     for i=1:r3
% % %         for j=1:c3
% % %             if I3(i,j,1)>0  %%在行的方向上找到第一个不为0 的点的坐标[rr2,cc2]
% % %                 rr2=i;          %在列方向上定位最左边 （左边缘）
% % %                 %         cc2=j;
% % %                 flag=1;
% % %                 break;
% % %             end
% % %         end
% % %         if flag==1
% % %             break;
% % %         end
% % %     end
% % %     r=rr2;
% % %     c=cc1;
% % %     
% % %     rstart=r;
% % %     %因为旋转后的图像是对称的，所以画布两端非图像部分的宽度是一样的
% % %     rend=r3-r+1;
% % %     cstart=c;
% % %     cend=c3-c+1;
% % %     I4=I3(rstart+1:rend-1,cstart+1:cend-1,:);
% % %     subplot(2,3,4),imshow(I4),title('裁剪边缘');
% % %     correctimg=imresize(I4,[512 512]);%%%对裁边的图像进行缩放到规则尺寸
% % %     psnr=colorpsnr(watermarkedimg,correctimg);%%计算校准图像与原始图像Host的相似度PSNR
% % %     disp(['测得的旋转角度：',num2str(angle),'    PSNR=',num2str(psnr)]);    %%求的旋转角度
% % %     subplot(2,3,5),imshow(correctimg);title(['规范尺寸 PSNR=',num2str(psnr,-4)]);