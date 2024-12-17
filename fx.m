clc;
 im={'lena.jpg'}
 blocksize=4;
 T=160;
      result(T,1)=T;
% T=0.012 
%ʵ����������block=4��4���к�ǿ�Ĳ��ɼ��ԣ�
%��Խ�����ԣ�block=8��8���нϺõ�³���ԣ�Ч����ǰ������޼���
%��ˮӡ�����ϴ�ʱ���Կ����ý�С�ķֿ�ߴ�õ�����ֿ���������ˮӡ�� 
%�ı�ֱ�ı�U��
%H='lena.jpg';%;'peppers.jpg';%'baboon.jpg';%'F16.jpg';%
H=im{1,1};
%H='lena.jpg';
W='3073232.tif';
%W='200100732.jpg';
%W='colorwatermark32.bmp';
%W='beijing1.bmp';
%W='ldu3232.jpg';
Host=imread(H);
Water=imread(W);
lenw=size(Water,1)*size(Water,2)*8;
% figure(1),subplot(121),imshow(Host),title('Host image');
% subplot(122),imshow(Water),title('watermark image');
%ˮӡArnold�任
%--------------------------------------------------
%ԭʼˮӡ���� �����Ӻ����γɷֲ�ˮӡlevelwatermark
%water307=imread('colorwatermark32.bmp');
tic;

for level=1:3
    levelwatermark(level,:)=gainlevelwatermark(Water(:,:,level));
end
%--------------------------------------------------
%ˮӡǶ��
psnrval=0;
 for level=1:3
% level=3
Hostlevel=double(Host(:,:,level));
Hostlevelback=Hostlevel;
nlevels =1;%[0, 1, 3] ;        % Decomposition level
pfilter = 'maxflat' ;%'9-7';%              % Pyramidal filter
dfilter = 'dmaxflat7' ; %  'pkva';           % Directional filter

% Nonsubsampled Contourlet decomposition
coeffs = nsctdec( Hostlevel, nlevels, dfilter, pfilter );

Lowpass=coeffs{1,1};%%save in su1

Waterold=levelwatermark(level,:);
[Hr,Hc]=size(Lowpass);
waterposition=1;

for  i=1:Hr/blocksize%������ͼ��ֿ飬���п��Էֳɼ�����
        
        if mod(i,2)==1 && mod(level,2)==1%���ݲ�������������һ������ѡ�鲼��
            star=1;
        else
            star=2;
        end
        
      for  j=star:2:Hc/blocksize%�ȶ��н��л��֣��а���4*4�ֿ�ԭ����Էֶ��ٴΣ����ѡ�У�������У�����Ϊ2
        watermark=Waterold(1,waterposition);%����ʱ��Ƕ���ˮӡ�ĸ�������ˮӡλ
        block=Lowpass((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);%�ֿ�λ��
        h=sum((sum(block))');%������ϵ����
        
        if (watermark=='1') %ˮӡǶ��ֵ����1��ʱ�������������¡�
           if mod(h,T)<0.25*T
               h1=h-mod(h,T)-0.25*T;
           else
               h1=h-mod(h,T)+0.75*T;
           end
        else %ˮӡǶ��ֵ����0��ʱ��
           if mod(h,T)>=0.75*T
               h1=h-mod(h,T)+1.25*T;
           else
                h1=h-mod(h,T)+0.25*T;
           end
        end
     d=(h1-h)/20; 
        
     for m=1:4
            for n=1:4
               if (m==2&&n==2)||(m==2&&n==3)||(m==3&&n==2)||(m==3&&n==3)
                   block(m,n)=block(m,n)+2*d;
               else
                   %if (m==1&&n==1)|| (m==1&&n==4)||(m==4&&n==1)||(m==4&&n==4)
                   block(m,n)=block(m,n)+d;
%                    end
               end
            end
     end
     Lowpass((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize)=block;%���ֿ�ֵ���¸���ԭ���ķֿ�
     waterposition=waterposition+1;%ˮӡǶ��λ�ã�Ƕ��������1
     if waterposition>lenw %%%==2799%%�ж��Ƿ����е�ˮӡ���ж�Ƕ���ȥ
            break;
     end
    end     %�����ѭ�����˽��� 
     if waterposition>lenw %%==5778==5578==2799
           break;
     end
end %��һ��ˮӡͼ��Ƕ�����
coeffs{1,1}=Lowpass;
Hostlevel= nsctrec(coeffs, dfilter, pfilter) ;
watermarkedim(:,:,level)=Hostlevel;
imwrite(uint8(Hostlevel),'www.bmp');
imwrite(Host(:,:,level),'hhh.bmp');
 end
imwrite(uint8(watermarkedim),'watermarked.bmp');
psnrval=colorpsnr(watermarkedim,Host);
ssimval=colorssim(watermarkedim,Host);
subplot(131),imshow(Host),title('Original image');
subplot(132),imshow('watermarked.bmp'),title([' PSNR=',num2str(psnrval),' SSIM=',num2str(ssimval),' T=',num2str(T)]);
result(T,2)=psnrval;
result(T,3)=ssimval;
toc;
 tic; 
%��ȡˮӡ
watermarkedim=double(imread('watermarked.bmp'));
for level=1:3
    wa=watermarkedim(:,:,level);
    Waterold=levelwatermark(level,:);
% Nonsubsampled Contourlet decomposition
    coeffs = nsctdec( wa, nlevels, dfilter, pfilter );
    wa=coeffs{1,1};%%su3
    waterposition=1;
     for  i=1:Hr/blocksize%������ͼ��ֿ飬���п��Էֳɼ�����
         if mod(i,2)==1 && mod(level,2)==1%���ݲ�������������һ������ѡ�鲼��
            star=1;
        else
            star=2;
        end
        
      for  j=star:2:Hc/blocksize%�ȶ��н��л��֣��а���4*4�ֿ�ԭ����Էֶ��ٴΣ����ѡ�У�������У�����Ϊ2
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);%�ֿ�λ��
        h2=sum((sum(block))');
        if  mod(h2,T)>0.5*T
            watermark='1';
        else
            watermark='0';
        end
        ExWater(1,waterposition)=watermark;%����ȡ�����Ķ�����ֵ������ȡ������ ExWater����
        waterposition=waterposition+1;
       if waterposition>lenw %==5578%==2799
           break;
       end
      end
     end
     wateronedim(level,:)=ExWater;
end

% %--------------------------------------------------
%ˮӡ�ָ�
extractlevelwatermark=watermarkrestore(wateronedim);
toc;
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=nc(uint8(extractlevelwatermark),uint8(Water));
subplot(133),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval)]);
name=[H(1:length(H)-4),'Unattack',W(1:length(W)-4),'.fig'];
%saveas(hand,name);
result(T,4)=ncval;

%%�㷨����
watermarkedim=imread('watermarked.bmp');
 for Q=10:10:100
   % Q=70;
    for level=1:3
    imwrite(watermarkedim(:,:,level),'temp.bmp','jpg','quality',Q);%����JPEGѹ������
    wb=imread('temp.bmp');
    ttttt(:,:,level)=wb;
    coeffs1 = nsctdec(wb, nlevels, dfilter, pfilter );
    Lowpass1=coeffs1{1,1};%%save in su1
    wa=double(Lowpass1);
    waterposition=1;
    for i=1:Hr/blocksize
        if mod(i,2)==1 && mod(level,2)==1
            star=1;
        else
            star=2;
        end
    for j=star:2:Hc/blocksize
        blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
         h2=sum((sum(blocknew))');
        if  mod(h2,T)>0.5*T
            watermark='1';
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
%    extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
 end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
hand=figure(Q/10+2);
subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark JPEG nc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenaJPEG', num2str(Q/10),'307','.fig'];
saveas(hand,name)
%result(ss,3+Q/10)=ncval;
 %ss2=3+Q/10;
 end
 %%
 %JPEG2000
%  watermarkedim=imread('watermarked.bmp');
%  for Q=1:10
%     for level=1:3
%     imwrite(watermarkedim(:,:,level),'temp.bmp','jp2','compressionratio',Q);
%     wb=imread('temp.bmp');
%     ttttt(:,:,level)=wb;
%     coeffs1 = nsctdec(wb, nlevels, dfilter, pfilter );
%     Lowpass1=coeffs1{1,1};%%save in su1
%     wa=double(Lowpass1);
%     waterposition=1;
%     for i=1:Hr/blocksize
%         if mod(i,2)==1 && mod(level,2)==1
%             star=1;
%         else
%             star=2;
%         end
%     for j=star:2:Hc/blocksize
%         blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         h3=sum((sum(blocknew))');
%         if  mod(h3,T)>0.5*T
%             watermark='1';
%         else
%             watermark='0';
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
% %     extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
% end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark JPEG2000 nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenaJPEG2000', num2str(Q),'307','.fig'];
% % saveas(hand,name)
% % 
% % result(ss,ss2+Q)=ncval;
%  end
 watermarkedim=imread('watermarked.bmp');
 for Q=1:5
    for level=1:3
             wb=imnoise(watermarkedim(:,:,level),'gaussian',0,0.001*Q);
    ttttt(:,:,level)=wb;
    coeffs1 = nsctdec(wb, nlevels, dfilter, pfilter );
    Lowpass1=coeffs1{1,1};%%save in su1
    wa=double(Lowpass1);
    waterposition=1;
    for i=1:Hr/blocksize
        if mod(i,2)==1 && mod(level,2)==1
            star=1;
        else
            star=2;
        end
    for j=star:2:Hc/blocksize
       blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
       h3=sum((sum(blocknew))');
        if  mod(h3,T)>0.5*T
            watermark='1';
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

for level=1:3
    extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
    %extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
hand=figure(Q+2);
subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark gaussiannc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenagaussian', num2str(Q),'by 0.1 right','.fig'];
% saveas(hand,name)
% result(ss,ss2+Q)=ncval;
 end
% watermarkedim=imread('watermarked.bmp');
%  for Q=2:5
%     for level=1:3
%         wb=imresize(watermarkedim(:,:,level),1/Q); 
%         wb=imresize(wb,[512 512]);   
%          ttttt(:,:,level)=wb;
%     coeffs1 = nsctdec(wb, nlevels, dfilter, pfilter );
%     Lowpass1=coeffs1{1,1};%%save in su1
%     wa=double(Lowpass1);
%     waterposition=1;
%        for i=1:Hr/blocksize
%        if mod(i,2)==1 && mod(level,2)==1
%             star=1;
%         else
%             star=2;
%         end
%     for j=star:2:Hc/blocksize
%         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%        blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         h3=sum((sum(blocknew))');
%         if  mod(h3,T)>0.5*T
%             watermark='1';
%         else
%             watermark='0';
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
 
%ˮӡ�ظ�
% 
% for level=1:3
%  extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
% %extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
% end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark Resize small nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenaSCaling','by 111', num2str(Q),'.fig'];
% % saveas(hand,name)
% % result(ss,ss2+Q)=ncval;
%  end
% watermarkedim=imread('watermarked.bmp');
%  for Q=1:5
%     for level=1:3
%        wb=imnoise(watermarkedim(:,:,level),'salt & pepper',0.01*Q);  %salt & pepper
%      ttttt(:,:,level)=wb;
%     coeffs1 = nsctdec(wb, nlevels, dfilter, pfilter );
%     Lowpass1=coeffs1{1,1};%%save in su1
%     wa=double(Lowpass1);
%     waterposition=1;
%     for i=1:Hr/blocksize
%        if mod(i,2)==1 && mod(level,2)==1
%             star=1;
%         else
%             star=2;
%         end
%     for j=star:2:Hc/blocksize
%        blocknew=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         h3=sum((sum(blocknew))');
%         if  mod(h3,T)>0.5*T
%             watermark='1';
%         else
%             watermark='0';
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
%     %extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
% end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark salt nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenaSalt', num2str(Q),'by0.02','.fig'];
% % saveas(hand,name)
% % result(ss,ss2+Q)=ncval;
%  end
