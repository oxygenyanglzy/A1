function rotateimg=RotatingAttack(img,angle)
%��ת����ͼ�񣬲�У����2019-11-24���
I2=imrotate(img,angle);
% M=512;N=512;
% if angle==0
%     rotateimg=img;
%     disp('û�ж�ͼ�������ת����!');
% else
%     [r2,c2,~]=size(I2);
%     %%��ת�ǶȲ���% ��Ƕ���ת����
%     figure(10);
%     subplot(2,3,1),imshow(uint8(img));title('lena');
%     for j=1:c2
%         for i=1:r2
%             if I2(i,j,1)>0 %������תͼ�������Ե�Ľ���
%                 ry=i;           %�ڴ�ֱ���������ȶ�λ
%                 rx=r2-ry;
%                 break;
%             end
%         end
%     end
%     angle=round(atan(rx/ry)*180/pi);%������ת�Ƕ�
%     
%      I3=imrotate(I2,-angle);%%%��������ǶȽ���������ת
%     subplot(2,3,2),imshow(uint8(I2));title(['rotated image ','angle=',num2str(angle)]);
%     subplot(2,3,3),imshow(uint8(I3));title('Inverse rotation image');
%     %%----------------�ü���Ե-----------------
%     %%-------�ü���Ե��ȥ���ڱ߲����ŵ�ָ����С[M,N]-------%%
% [f,c]=find(I3(:,:,3));
% minf=min(f);maxf=max(f);
% minc=min(c);maxc=max(c);
% 
% %����ͼƬ��߷���Ԫ�صĸ���
% count=0;
% for j=minc:maxc
%     for i=minf:maxf
%         if I3(i,j,3)>=30 && I3(i,j,2)>=30 && I3(i,j,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(maxf-minf)/2  %�������Ԫ�ذ���������<30��Ԫ�أ���˵��������Ч������+1��ֱ����ЧΪֹ
%         minc=minc+1;
%     else
%         break;
%     end
% end
% %����ͼƬ�ұ߷���Ԫ�صĸ���
% count=0;
% for j=minc:maxc
%     jj=maxc+minc-j;
%     for i=minf:maxf
%         if I3(i,jj,3)>=30 && I3(i,jj,2)>=30 && I3(i,jj,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(maxf-minf)/2 %���β��Ԫ�ذ���������<30��Ԫ�أ���˵��������Ч������-1��ֱ����ЧΪֹ
%         maxc=maxc-1;
%     else
%         break;
%     end
% end
% %����ͼƬ�ϱ߷���Ԫ�صĸ���
% count=0;
% for i=minf:maxf
%     for j=minc:maxc
%         if I3(i,j,3)>=30 && I3(i,j,2)>=30 && I3(i,j,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(maxc-minc)/2 %�������Ԫ�ذ���������<30��Ԫ�أ���˵��������Ч������+1��ֱ����ЧΪֹ
%         minf=minf+1;
%     else
%         break;
%     end
% end
% %����ͼƬ�±߷���Ԫ�صĸ���
% count=0;
% for i=minf:maxf
%     ii=minf+maxf-i;
%     for j=minc:maxc
%         if I3(ii,j,3)>=30 && I3(ii,j,2)>=30 && I3(ii,j,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(maxc-minc)/2 %���β��Ԫ�ذ���������<30��Ԫ�أ���˵��������Ч������-1��ֱ����ЧΪֹ
%         maxf=maxf-1;
%     else
%         break;
%     end
% end
% I4=I3(minf:maxf,minc:maxc,:);
% correctimg=imresize(I4,[M N]);%%%�Բñߵ�ͼ��������ŵ�����ߴ�
% subplot(2,3,4),imshow(uint8(I4)),title('�ü���Ե');
% correctimg=imresize(I4,[512 512]);%%%�Բñߵ�ͼ��������ŵ�����ߴ�
% psnr=colorpsnr(watermarkedimg,correctimg);%%����У׼ͼ����ԭʼͼ��Host�����ƶ�PSNR
% disp(['��õ���ת�Ƕȣ�',num2str(angle),'    PSNR=',num2str(psnr)]);    %%�����ת�Ƕ�
% subplot(2,3,5),imshow(uint8(correctimg));title(['�淶�ߴ� PSNR=',num2str(psnr,-4)]);
% end
rotateimg=I2;
% figure;
% imshow(uint8(rotateimg));
% imwrite(uint8(I2),['Lena-Rotation',num2str(angle),'.jpg']);
end
% % % [r3,c3,~]=size(I3);
% % %     flag=0; %��λͼ��������
% % %     for j=1:c3
% % %         for i=1:r3
% % %             if I3(i,j,1)>0  %%���еķ������ҵ���һ����Ϊ0 �ĵ������[rr1,cc1]
% % %                 %         rr1=i;          %���з����϶�λ�����
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
% % %     flag=0;%��λͼ��������
% % %     for i=1:r3
% % %         for j=1:c3
% % %             if I3(i,j,1)>0  %%���еķ������ҵ���һ����Ϊ0 �ĵ������[rr2,cc2]
% % %                 rr2=i;          %���з����϶�λ����� �����Ե��
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
% % %     %��Ϊ��ת���ͼ���ǶԳƵģ����Ի������˷�ͼ�񲿷ֵĿ����һ����
% % %     rend=r3-r+1;
% % %     cstart=c;
% % %     cend=c3-c+1;
% % %     I4=I3(rstart+1:rend-1,cstart+1:cend-1,:);
% % %     subplot(2,3,4),imshow(I4),title('�ü���Ե');
% % %     correctimg=imresize(I4,[512 512]);%%%�Բñߵ�ͼ��������ŵ�����ߴ�
% % %     psnr=colorpsnr(watermarkedimg,correctimg);%%����У׼ͼ����ԭʼͼ��Host�����ƶ�PSNR
% % %     disp(['��õ���ת�Ƕȣ�',num2str(angle),'    PSNR=',num2str(psnr)]);    %%�����ת�Ƕ�
% % %     subplot(2,3,5),imshow(correctimg);title(['�淶�ߴ� PSNR=',num2str(psnr,-4)]);