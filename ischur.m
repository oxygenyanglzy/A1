clear all;
blocksize=4;
sgm=25;
% %�����������ң����϶��µ�˳�������ѡ��Ƕ��� ����������DC����������
H='F16.jpg';
%W='ldu3232.jpg';
W='3073232.tif';
%W='200100732.jpg';
Host=imread(H);%��������ͼƬ
Water=imread(W);%����ˮӡ
lenw=size(Water,1)*size(Water,2)*8;%ˮӡ��������*������32*32*8ˮӡͼ����ռ��bitλ��
figure(1),subplot(121),imshow(Host),title('Host image');
subplot(122),imshow(Water),title('watermark image');%��ʾˮӡͼƬ������ͼƬ
%--------------------------------------------------
%ԭʼˮӡ���� �����Ӻ����γɷֲ�ˮӡlevelwatermark
for level=1:3
    temp1=Water(:,:,level);%����ֲ���ȡ
    afterarnold=Arnold(temp1,6,0);%��ˮӡͼ��ֲ����ң�������6�Σ�����
    levelwatermark(level,:)=gainlevelwatermark(afterarnold);%��ü��ܺ�ķֲ�ͼ����󣬲�ת���ɶ�����������ˮӡ��Ϣ
end
%--------------------------------------------------
%ˮӡǶ��
tic;
psnrval=0;
 for level=1:3
Hostlevel=double(Host(:,:,level));%�Ƚ�������ͼ�񣩷ֲ㲢ת����double���͵Ķ�ά����
% Hostlevelback=Hostlevel;
Waterold=levelwatermark(level,:);%����Ӧ��ˮӡlevel���������и���waterold��ˮӡ�������㣬һ��һ���Ƕ��
[Hr,Hc]=size(Hostlevel);%����Hostlevel��ά�������������
waterposition=1;%Ƕ��ˮӡλ�ã���ʼֵΪ1
for i=1:Hr/blocksize%������ͼ��ֿ飬���п��Էֳɼ�����
        if mod(i,2)==1 && mod(level,2)==1%���ݲ�������������һ������ѡ�鲼��
            star=1;
        else
            star=2;
        end
    for j=star:2:Hc/blocksize%�ȶ��н��л��֣��а���4*4�ֿ�ԭ����Էֶ��ٴΣ����ѡ�У�������У�����Ϊ2
        watermark=Waterold(1,waterposition);%����ʱ��Ƕ���ˮӡ����ֵ����ˮӡλ
        block=Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);%�ֿ�λ��
        [u,s]=schur(block);
        if level==2
                sgm1=sgm*0.94;
        else
            if level==3
                sgm1=sgm*0.78;
            else
                sgm1=sgm;
            end
        end
         lmt=floor(s(1,1)/sgm1);
        if (watermark=='1') %ˮӡǶ��ֵ����1��ʱ�������������¡�
           if mod(lmt+1,2)==1
                     s1=(lmt+1.5)*sgm1;
           else
                      s1=(lmt+0.5)*sgm1;
            end
        else %ˮӡǶ��ֵ����0��ʱ��
           if mod(lmt,2)==1
                       s1=(lmt+1.5)*sgm1;
           else
                       s1=(lmt+0.5)*sgm1;
           end
        end
        s(1,1)=s1;
        blocknew=u*s*u';
        if max(max(blocknew))>255%�жϿ�����ֵ�Ƿ�Խ��
           block=Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
            if (watermark=='1') %ˮӡǶ��ֵ����1��ʱ�������������¡�
               if mod(lmt+1,2)==1
                     s1=(lmt-0.5)*sgm1;
               else
                     s1=(lmt-1.5)*sgm1;
               end
           else %ˮӡǶ��ֵ����0��ʱ��
               if mod(lmt,2)==1
                       s1=(lmt-0.5)*sgm1;
               else
                       s1=(lmt-1.5)*sgm1;
               end
           end
        s(1,1)=s1;
        blocknew=u*s*u';
       end
        Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize)=blocknew;%���ֿ�ֵ���¸���ԭ���ķֿ�
        waterposition=waterposition+1;%ˮӡǶ��λ�ã�Ƕ��������1
        if waterposition>lenw %%%==2799%%�ж��Ƿ����е�ˮӡ���ж�Ƕ���ȥ
            break;
        end
    end
     if waterposition>lenw %%==5778==5578==2799
           break;
     end
    end%��һ��ˮӡͼ��Ƕ�����
watermarkedim(:,:,level)=Hostlevel;%����һ�������ͼ���ά���󸳸�
imwrite(uint8(Hostlevel),'www.bmp');%��Hostlevel��Ƕ��ˮӡ�������ͼ�񣩣�ԭ��Ϊdouble���͵ģ�ת����uint8,����ͼƬ.png��ʽ
imwrite(Host(:,:,level),'hhh.bmp');%����������һ���ͼ��
psnrval=psnrval+PSNR('www.bmp','hhh.bmp');%���������һ�飬��������ͼƬ֮������ƶ�
 end%RGB����ͼ�����Ƕ�����
hfx=watermarkedim;
imwrite(uint8(watermarkedim),'watermarked.bmp');%����Ƕ��ˮӡ��ͼƬΪwatermarked.bmp�����ڹ�������ͼƬ��
%psnrval=PSNR(H,'watermarked.bmp');
hand=figure(2),subplot(131),imshow(Host),title('Original image');
subplot(132),imshow('watermarked.bmp'),title([' unattack PSNR=',num2str(psnrval/3)]);%��ʾǶ����ɺ��ˮӡͼƬ
toc;
%--------------------------------------------------
%ˮӡ��ȡ
for level=1:3
    wa=watermarkedim(:,:,level);
    waterposition=1;
    for i=1:Hr/blocksize
        if mod(i,2)==1 && mod(level,2)==1
            star=1;
        else
            star=2;
        end
       for j=star:2:Hc/blocksize
        blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u1,s1]=schur(blocknew);
        if level==2
                sgm2=sgm*0.94;
        else
            if level==3
                sgm2=sgm*0.78;
             else
                sgm2=sgm;
            end
        end
        lmt1=floor(s1(1,1)/sgm2); 
        if  mod(lmt1,2)==1
               watermark= '1';
        else
               watermark='0';
        end
        
        ExWater(1,waterposition)=watermark;
        waterposition=waterposition+1;
        if waterposition>lenw
            break;
        end
     end

     if waterposition>lenw %%==5578
           break;
     end
    end
     wateronedim(level,:)=ExWater;%����Ӧ��ˮӡֵ������Ӧ��R,G,B�㣬
end
% %--------------------------------------------------
%ˮӡ�ָ�

for level=1:3
    extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%��ˮӡ����������ת������Ӧ������ֵ��
    extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);%�����Һ������ֵ����Arnold�任ת�����
    
 end
% wateronedim=reshape(ExWater,size(Water,1),size(Water,2));

imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=nc(uint8(extractlevelwatermark),uint8(Water));%�Ƚ�ˮӡͼ��ncֵ
subplot(133),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval)]);

%%�㷨����
%%
 watermarkedim=imread('watermarked.bmp');
 for Q=10:10:100
    for level=1:3
    imwrite(watermarkedim(:,:,level),'temp.bmp','jpg','quality',Q);%����JPEGѹ������
    wb=imread('temp.bmp');
    ttttt(:,:,level)=wb;
    wa=double(wb);
    waterposition=1;
    for i=1:Hr/blocksize
        if mod(i,2)==1 && mod(level,2)==1
            star=1;
        else
            star=2;
        end
    for j=star:2:Hc/blocksize
        blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u1,s1]=schur(blocknew);
        if level==2
                sgm2=sgm*0.94;
        else
            if level==3
                sgm2=sgm*0.78;
             else
                sgm2=sgm;
            end
        end
        lmt1=floor(s1(1,1)/sgm2); 
        if  mod(lmt1,2)==1
               watermark= '1';
        else
               watermark='0';
        end
        
        ExWater(1,waterposition)=watermark;
        waterposition=waterposition+1;
        if waterposition>lenw
            break;
        end
              
    end
     if waterposition>lenw
           break;
     end
    end
    wateronedim(level,:)=ExWater;
 end
 
%%ˮӡ�ָ�

for level=1:3
   extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
    extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
 end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
hand=figure(Q/10+2);
subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark JPEG nc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenaJPEG', num2str(Q/10),'307','.fig'];
% saveas(hand,name)
% result(ss,3+Q/10)=ncval;
end
%%JPEG2000
%  watermarkedim=imread('watermarked.bmp');
%  for Q=1:10
%     for level=1:3
%         imwrite(watermarkedim(:,:,level),'temp.bmp','jp2','compressionratio',Q);
%     wb=imread('temp.bmp');
%     ttttt(:,:,level)=wb;
%     wa=double(wb);
%     waterposition=1;
%     for i=1:Hr/blocksize
%         if mod(i,2)==1 && mod(level,2)==1
%             star=1;
%         else
%             star=2;
%         end
%     for j=star:2:Hc/blocksize
%        blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         [u1,s1]=schur(blocknew);
%         if level==2
%                 sgm2=sgm*0.94;
%         else
%             if level==3
%                 sgm2=sgm*0.78;
%              else
%                 sgm2=sgm;
%             end
%         end
%           lmt1=floor(s1(1,1)/sgm2); 
%         if  mod(lmt1,2)==1
%                watermark= '1';
%         else
%                watermark='0';
%         end
%         ExWater(1,waterposition)=watermark;
%         waterposition=waterposition+1;
%         if waterposition>lenw
%             break;
%         end
%               
%     end
%      if waterposition>lenw
%            break;
%      end
%     end
%     wateronedim(level,:)=ExWater;
%      
% end
%  
% %ˮӡ�ظ�
% 
% for level=1:3
%      extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
%     extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
% end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark JPEG2000 nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenaJPEG2000', num2str(Q),'307','.fig'];
%  end
%%��������
% watermarkedim=imread('watermarked.bmp');
%  for Q=1:5
%     for level=1:3
%              wb=imnoise(watermarkedim(:,:,level),'gaussian',0,0.001*Q);
%     ttttt(:,:,level)=wb;
%     wa=double(wb);
%     waterposition=1;
%     for i=1:Hr/blocksize
%         if mod(i,2)==1 && mod(level,2)==1
%             star=1;
%         else
%             star=2;
%         end
%     for j=star:2:Hc/blocksize
%        blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         [u1,s1]=schur(blocknew);
%         if level==2
%                 sgm2=sgm*0.94;
%         else
%             if level==3
%                 sgm2=sgm*0.78;
%              else
%                 sgm2=sgm;
%             end
%         end
%         lmt1=floor(s1(1,1)/sgm2); 
%         if  mod(lmt1,2)==1
%                watermark= '1';
%         else
%                watermark='0';
%         end
%         ExWater(1,waterposition)=watermark;
%         waterposition=waterposition+1;
%         if waterposition>lenw
%             break;
%         end
%               
%     end
%      if waterposition>lenw
%            break;
%      end
%     end
%     wateronedim(level,:)=ExWater;
%      
% end
 
% ˮӡ�ظ�
% 
% for level=1:3
%     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
%     extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
% end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark gaussiannc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenagaussian', num2str(Q),'by 0.1 right','.fig'];
% end
%%��������
% watermarkedim=imread('watermarked.bmp');
%  for Q=1:5
%     for level=1:3
%        wb=imnoise(watermarkedim(:,:,level),'salt & pepper',0.02*Q);  %salt & pepper
%        ttttt(:,:,level)=wb;
%     wa=double(wb);
%     waterposition=1;
%     for i=1:Hr/blocksize
%        if mod(i,2)==1 && mod(level,2)==1
%             star=1;
%         else
%             star=2;
%         end
%     for j=star:2:Hc/blocksize
%         blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         [u1,s1]=schur(blocknew);
%         if level==2
%                 sgm2=sgm*0.94;
%         else
%             if level==3
%                 sgm2=sgm*0.78;
%              else
%                 sgm2=sgm;
%             end
%         end
%         lmt1=floor(s1(1,1)/sgm2); 
%         if  mod(lmt1,2)==1
%                watermark= '1';
%         else
%                watermark='0';
%         end
%         ExWater(1,waterposition)=watermark;
%         waterposition=waterposition+1;
%         if waterposition>lenw
%             break;
%         end
%      end
%      if waterposition>lenw
%            break;
%      end
%     end
%     wateronedim(level,:)=ExWater;
%  end
%  
% %ˮӡ�ظ�
% 
% for level=1:3
%     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
%     extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
% end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark salt nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenaSalt', num2str(Q),'by0.02','.fig'];
%  end
%%��ֵ�˲�
% watermarkedim=imread('watermarked.bmp');
%  for Q=1:5
%     for level=1:3
%       wb=medfilt2(watermarkedim(:,:,level),[Q,Q]);
%       ttttt(:,:,level)=wb;
%     wa=double(wb);
%     waterposition=1;
%     for i=1:Hr/blocksize
% 
%         if mod(i,2)==1 && mod(level,2)==1
%             star=1;
%         else
%             star=2;
%         end
%     for j=star:2:Hc/blocksize
%       blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         [u1,s1]=schur(blocknew);
%         if level==2
%                 sgm2=sgm*0.94;
%         else
%             if level==3
%                 sgm2=sgm*0.78;
%              else
%                 sgm2=sgm;
%             end
%         end
%         lmt1=floor(s1(1,1)/sgm2); 
%         if  mod(lmt1,2)==1
%                watermark= '1';
%         else
%                watermark='0';
%         end
%         ExWater(1,waterposition)=watermark;
%         waterposition=waterposition+1;
%         if waterposition>lenw
%             break;
%         end
%               
%     end
%      if waterposition>lenw
%            break;
%      end
%     end
%     wateronedim(level,:)=ExWater;
%      
% end
%  
% %ˮӡ�ظ�
% 
% for level=1:3
%       extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
%     extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);;
% end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark Median nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenaMedian', num2str(Q),'by','.fig'];
%  end
%%��ͨ�˲�
% watermarkedim=imread('watermarked.bmp');
%  for Q=1:5
%     for level=1:3
%    Ib=watermarkedim(:,:,level);
% f=double(Ib);
% g=fft2(f);
% g=fftshift(g);
% [N1,N2]=size(g);
% d0=100;
% n1=fix(N1/2);
% n2=fix(N2/2);
% for i=1:N1
%     for j=1:N2
%         d=sqrt((i-n1)^2+(j-n2)^2);
%         h=1/(1+0.414*(d/d0)^(2*Q));
%         resul(i,j)=h*g(i,j);
%          end
% end
% resul=ifftshift(resul);
% X2=ifft2(resul);
% X3=uint8(real(X2));
%     ttttt(:,:,level)=X3;
%     wa=double(X3);
%     waterposition=1;
%     for i=1:Hr/blocksize
%        if mod(i,2)==1 && mod(level,2)==1
%             star=1;
%         else
%             star=2;
%         end
%     for j=star:2:Hc/blocksize
%           blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         [u1,s1]=schur(blocknew);
%         if level==2
%                 sgm2=sgm*0.94;
%         else
%             if level==3
%                 sgm2=sgm*0.78;
%              else
%                 sgm2=sgm;
%             end
%         end
%         lmt1=floor(s1(1,1)/sgm2); 
%         if  mod(lmt1,2)==1
%                watermark= '1';
%         else
%                watermark='0';
%         end
%         ExWater(1,waterposition)=watermark;
%         waterposition=waterposition+1;
%         if waterposition>lenw
%             break;
%         end
%        end
%      if waterposition>lenw
%            break;
%      end
%     end
%     wateronedim(level,:)=ExWater;
%     end
%  %ˮӡ�ظ�
% for level=1:3
%  extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
%    extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
% end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark Butterlow-pass nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenaButterLowpass', 'by',num2str(Q),'.fig'];
% end
%%resize����
% watermarkedim=imread('watermarked.bmp');
%  for Q=2:5
%     for level=1:3
%         wb=imresize(watermarkedim(:,:,level),1/Q); 
%         wb=imresize(wb,[512 512]);   
%         ttttt(:,:,level)=wb;
%         wa=double(wb);
%         waterposition=1;
%        for i=1:Hr/blocksize
%        if mod(i,2)==1 && mod(level,2)==1
%             star=1;
%         else
%             star=2;
%         end
%     for j=star:2:Hc/blocksize
%          blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         [u1,s1]=schur(blocknew);
%         if level==2
%                 sgm2=sgm*0.94;
%         else
%             if level==3
%                 sgm2=sgm*0.78;
%              else
%                 sgm2=sgm;
%             end
%         end
%         lmt1=floor(s1(1,1)/sgm2); 
%         if  mod(lmt1,2)==1
%                watermark= '1';
%         else
%                watermark='0';
%         end
%         ExWater(1,waterposition)=watermark;
%         waterposition=waterposition+1;
%         if waterposition>lenw
%             break;
%         end
%               
%     end
%      if waterposition>lenw
%            break;
%      end
%     end
%     wateronedim(level,:)=ExWater;
%      
% end
%  
% %ˮӡ�ظ�
% 
% for level=1:3
%  extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
% extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
% end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark Resize small nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenaSCaling','by 111', num2str(Q),'.fig'];
%  end
%%���й���
watermarkedim=imread('watermarked.bmp');
wb=watermarkedim;
 for Q=1:4
%    wb(1:64*Q,:,:)=255;
   wb(1:64*Q,1:64*Q,:)=255;
   for level=1:3
    ttttt(:,:,level)=wb(:,:,level);
    wa=double(wb(:,:,level));
    waterposition=1;
    for i=1:Hr/blocksize
       if mod(i,2)==1 && mod(level,2)==1
            star=1;
        else
            star=2;
        end
    for j=star:2:Hc/blocksize
        blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u1,s1]=schur(blocknew);
        if level==2
                sgm2=sgm*0.94;
        else
            if level==3
                sgm2=sgm*0.78;
             else
                sgm2=sgm;
            end
        end
        lmt1=floor(s1(1,1)/sgm2); 
        if  mod(lmt1,2)==1
               watermark= '1';
        else
               watermark='0';
        end
        ExWater(1,waterposition)=watermark;
        waterposition=waterposition+1;
        if waterposition>lenw
            break;
        end
              
    end
     if waterposition>lenw
           break;
     end
    end
    wateronedim(level,:)=ExWater;
     
end
 
% % ˮӡ�ظ�

for level=1:3
    tt=gainlevelimage(wateronedim(level,:));
    te=Arnold(tt,6,1);
    extimage(:,:,level)=te;
%     imshow(te);
 end
imwrite(ttttt,'temp1.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8( extimage),'extrwater1.bmp')
ncval=colornc(uint8( extimage),uint8(Water));
hand=figure(Q+2);
subplot(121),imshow( uint8( extimage)),title(['Extracted Watermark Crop nc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenaCrop64by',num2str(Q), '.fig'];
 end
 
