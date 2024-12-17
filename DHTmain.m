 clc
clear all;
blocksize=4;
load position;
R = 5;  % Maximum radius
N = 4;  % Number of sampling points
n = 0;
im={'baboon.jpg'};
k=1; 
for Qfac=745
  m=1;
H='lgthouse.bmp';
W='QQ.jpg';
Host=imread(H);
Water=imread(W);
lenw=size(Water,1)*size(Water,2)*8;
figure(1),subplot(121),imshow(Host),title('Host image');
subplot(122),imshow(Water),title('watermark image');
%水印Arnold变换
%--------------------------------------------------
%原始水印处理 调用子函数形成分层水印levelwatermark
%water307=imread('colorwatermark32.bmp');
%tic;
for level=1:3
    temp1=Water(:,:,level);%矩阵分层提取
    AfterHundun= Hundun(temp1,1);
    levelwatermark(level,:)=gainlevelwatermark(AfterHundun);
    % imshow(Hundun(temp1,1));
end
%--------------------------------------------------
%水印嵌入
%psnrval(1,k)=0;
        t1=0;
        t2=0;
        t1=cputime;
  for level=1:3
% level=3
Hostlevel=double(Host(:,:,level));
Hostlevelback=Hostlevel;  
fprintf("2123");
nlevels =1;%[0, 1, 3] ;        % Decomposition level
pfilter = 'maxflat' ;%'9-7';%              % Pyramidal filter
dfilter = 'dmaxflat7' ; %  'pkva';           % Directional filter
% Nonsubsampled Contourlet decomposition
% coeffs = nsctdec( Hostlevel, nlevels, dfilter, pfilter );
% Hostlevel=coeffs{1,1};%%save in su1

Waterold=levelwatermark(level,:);
[Hr,Hc]=size(Hostlevel);
waterposition=1;

 A1=0;
 A2=0;
 A3=0;
 A4=0;
for counter=1:lenw
    i=position(1,counter);%根据position取出应该的嵌入位置
    j=position(2,counter);
    watermark=Waterold(1,counter);%%根据水印序列取出水印位
    block=Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
    %将经过变换后的 v b
    h5 = @(r) block;
    [s,k1,r,I1,K1,R1,h2]=dht(h5,R,N,n);
    P = hadamard(4);
    P_inv = inv(P);
  s = s*P;
    lmt= round(s(1,1)/Qfac);
     if watermark=='1'
         wat=1;
     else
         wat=0;
     end
    Econd=xor(mod(lmt,2),wat);
    if Econd==0
        s(1,1)=lmt*Qfac+(Qfac/2);
    else
        s(1,1)=lmt*Qfac-(Qfac/2);
    end
   s2=s*P_inv
        blocknew=idht(s2, I1, K1, R1);
     Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize)=blocknew; 
end

% Reconstruct image
    watermarkedim(:,:,level)=Hostlevel;
    imwrite(uint8(Hostlevel),'www.bmp');
    imwrite(Host(:,:,level),'hhh.bmp');
    ssim(level)=colorssim(watermarkedim(:,:,level),Host(:,:,level));
  end
 %toc;
 time1=toc;
% disp(['q代码运行时间：', num2str(time1), ' 秒']);
 %toc;
   result(k,1)=cputime-t1;
imwrite(uint8(watermarkedim),'watermarked.bmp');
psnrval(k,1)=colorpsnr(watermarkedim,Host);
ssimval(k,1)=(ssim(1)+ssim(2)+ssim(3))/3;
hand=figure(2);
subplot(131),imshow(Host),title('Original image');
subplot(132),imshow('watermarked.bmp'),title([' PSNR=',num2str(psnrval(k,1)),' SSIM=',num2str(ssimval(k,1)),' Qfac=',num2str(Qfac)]);


%--------------------------------------------------
%水印提取
tic;
t2=cputime;
watermarkedim=double(imread('watermarked.bmp'));
for level=1:3
    wa=watermarkedim(:,:,level);
    Waterold=levelwatermark(level,:);
% Nonsubsampled Contourlet decomposition
    for counter=1:lenw
        i=position(1,counter);
        j=position(2,counter);
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
             h5 = @(r) block;
    [s1,k2,r1,I2,K2,R2,h3]=dht(h5,R,N,n);
      s1 = s1*P;
        Rext=fix(s1(1, 1)/Qfac);
        ext=mod(Rext,2);
        % end
if all(ext == 0)
    watermark = '0';
else
    watermark = num2str(mode(ext)); % 获取出现次数最多的值并转换为字符串
end
        if ext==0
            watermark='0';
        else
            watermark='1';
        end
        ExWater(1,counter)=watermark;
    end
    wateronedim(level,:)=ExWater;
end
time2=toc;
disp(['提取代码运行时间：', num2str(time2), ' 秒'])
 result(k,2)=cputime-t2;
%--------------------------------------------------
%水印恢复
for level=1:3
    extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
    extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来

 end
%extractlevelwatermark=watermarkrestore(wateronedim);
%toc;
time2=toc;
% disp(['提取代码运行时间：', num2str(time2), ' 秒']);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp');
ncval(k,1)=nc(uint8(extractlevelwatermark),uint8(Water));
subplot(133),imshow( uint8(extractlevelwatermark)),title(['111Extracted Watermark nc=',num2str(ncval(k,1))]);

ncval(k,1)=nc(uint8(extractlevelwatermark),uint8(Water));
subplot(133),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,1))]);
k=k+1;
  end
  

% % % 
% % % % % % % %% JPEG2000 测试 一  图像整体攻击提取效果较差
%  watermarkedim=imread('watermarked.bmp');
%  for QA=7
%   for level=1:3
%       imwrite(watermarkedim(:,:,level),'temp.bmp','jp2','compressionratio',QA);
%       wb=double(imread('temp.bmp'));
%       ttttt(:,:,level)=wb;
%       wa = wb;
%     for  counter=1:lenw
%         i=position(1,counter);
%         j=position(2,counter);
%         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%             h5 = @(r) block;
%      [s2,k2,r1,I2,K2,R2,h3]=dht(h5,R,N,n);
%      s2=s2*P;
%         Rext=fix(s2(1,1)/Qfac);
%         ext=mod(Rext,2);
%         if ext==0
%             watermark='0';
%         else
%             watermark='1';
%         end
%         ExWater(1,counter)=watermark;    
% 
%      end
%      wateronedim(level,:)=ExWater;
% 
% end
% 
%  total_bits = lenw; % 每个层次的总位数
% error_count = 0; % 错误位计数
% 
% for level = 1:3
%     for counter = 1:total_bits
%         original_bit = levelwatermark(level, counter);
%         extracted_bit = wateronedim(level, counter);
% 
%         if original_bit ~= extracted_bit
%             error_count = error_count + 1;
%         end
%     end
% end
% 
% ber = error_count / (total_bits * 3); % 总误码率
% disp(['误码率 (BER): ', num2str(ber)]);
% % % % % 水印回复
% for level=1:3
%     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
%     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% 
% end
% m=m+1;
% imwrite(ttttt,'temp.bmp');
% psnrval(k,m)=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(QA+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(QA),'  PSNR=',num2str(psnrval(k,m))]);
% name=[H(1:length(H)-4),'JPEG2000_', num2str(QA),'-',W(1:length(W)-4),'.fig'];
% saveas(hand,name);
% % result(1,QA)=psnrval;
% % result(2,QA)=ncval;
% % 
% end

% % 
% % 
% %  %  
% % % % 
% % % % %---------------------------------------------------------------------
%JPEG 测试 二  图像分层攻击提取效果较好
% % 
%  watermarkedim=imread('watermarked.bmp');
%  for QA=50%10:10:100
%     for level=1:3
%       imwrite(watermarkedim(:,:,level),'temp.bmp','jpg','quality',QA);
%       wb=double(imread('temp.bmp'));
%       ttttt(:,:,level)=wb;
% 
% % % Nonsubsampled Contourlet decomposition
%      coeffs = nsctdec( wb, nlevels, dfilter, pfilter );
%      wa=coeffs{1,1};%%su3
%      wa = wb;
%     for  counter=1:lenw
%          i=position(1,counter);
%          j=position(2,counter);
%          block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%             h5 = @(r) block;
%     [s2,k2,r1,I2,K2,R2,h3]=dht(h5,R,N,n);
%     s2=s2*P;
%         Rext=fix(s2(1,1)/Qfac);
%         ext=mod(Rext,2);
%         if ext==0
%             watermark='0';
%         else
%             watermark='1';
%         end
%         ExWater(1,counter)=watermark;        
%     end
%    wateronedim(level,:)=ExWater;
% 
%     end
% 
%  total_bits = lenw; % 每个层次的总位数
% error_count = 0; % 错误位计数
% 
% for level = 1:3
%     for counter = 1:total_bits
%         original_bit = levelwatermark(level, counter);
%         extracted_bit = wateronedim(level, counter);
% 
%         if original_bit ~= extracted_bit
%             error_count = error_count + 1;
%         end
%     end
% end
% 
% ber = error_count / (total_bits * 3); % 总误码率
% disp(['误码率 (BER): ', num2str(ber)]);
% % % % 水印回复
% 
% for level=1:3
%     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
%     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% 
%  end
% 
% imwrite(ttttt,'temp.bmp');
% psnrval(k,1)=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval(k,1)=nc(uint8(extractlevelwatermark),uint8(Water)); 
% % hand=figure(QA/10+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,1))]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(QA),'  PSNRr=',num2str(psnrval(k,1))]);
% name=[H(1:length(H)-4),'JPEG_', num2str(QA),'-',W(1:length(W)-4),'.fig'];
% saveas(hand,name);
% k=k+1;
%  end

% % % 高斯噪声
% a=imread('watermarked.bmp');
% for QQ=1
%     wb=imnoise(a,'salt & pepper',0.001*QQ);  %salt & pepper
%    % wb=imnoise(a,'gaussian',0,0.001*(QQ)); %gausian
%    imwrite(wb,'temp.bmp');
%  %  imshow(wb);
%   for level=1:3
%      ws=double(wb(:,:,level));
%      % coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
%      % wa=coeffs{1,1};%%su3
%      wa=ws;
% for  counter=1:lenw
%     i=position(1,counter);
%     j=position(2,counter);
% 
%         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         h5 = @(r) block;
%   [s2,k2,r1,I2,K2,R2,h3]=dht(h5,R,N,n);
%   s2=s2*P;
%         Rext=fix(s2(1,1)/Qfac);
%         ext=mod(Rext,2);
%         if ext==0
%             watermark='0';
%         else
%             watermark='1';
%         end
%         ExWater(1,counter)=watermark;  
%      end
%      wateronedim(level,:)=ExWater;
% 
% end
% % % 
% %水印回复
% 
%  total_bits = lenw; % 每个层次的总位数
% error_count = 0; % 错误位计数
% 
% for level = 1:3
%     for counter = 1:total_bits
%         original_bit = levelwatermark(level, counter);
%         extracted_bit = wateronedim(level, counter);
% 
%         if original_bit ~= extracted_bit
%             error_count = error_count + 1;
%         end
%     end
% end
% 
% ber = error_count / (total_bits * 3); % 总误码率
% disp(['误码率 (BER): ', num2str(ber)]);
% %水印回复
% for level=1:3
%     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
%     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% end
% psnrval(k,m)=colorpsnr(wb,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(QQ);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
% subplot(122),imshow(uint8(wb)),title([num2str(QQ),' PSNR=',num2str(psnrval(k,m))]);
% 
% name=[H(1:length(H)-4),'gaussian-', num2str(0.001*QQ),'-',W(1:length(W)-4),'.fig'];
% saveas(hand,name);
% m=m+1;
% end

%%%%%%%%%%%%%%%%%
% % watermarkedim=imread('watermarked.bmp');
%  for QQ=3:3
% 
%   %h = fspecial('average');%('gaussian');%, 50, 45);%%('average')%,'disk','gaussian',
% 
%   for level=1:3
%       %%wb = imfilter(watermarkedim(:,:,level), h);
%       wb=medfilt2(watermarkedim(:,:,level),[QQ,1]);
%        ttttt(:,:,level)=wb;
%     ws=double(wb);
%     waterposition=1;
%      % coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
%      % wa=coeffs{1,1};%%su3
%      wa=ws;
% for  counter=1:lenw
%     i=position(1,counter);
%     j=position(2,counter);
% 
%         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         h5 = @(r) block;
%     [s2,k2,r1,I2,K2,R2,h3]=dht(h5,R,N,n);
%         Rext=fix(s2(1,1)/Qfac);
%         ext=mod(Rext,2);
%         if ext==0
%             watermark='0';
%         else
%             watermark='1';
%         end
%         ExWater(1,counter)=watermark;  
% end
%      wateronedim(level,:)=ExWater;
% 
% end
% 
% %%%水印回复
% 
% extractlevelwatermark=watermarkrestore(wateronedim);
% 
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=nc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(QQ+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(QQ),'  PSNRr=',num2str(psnrval)]);
% name=[H(1:length(H)-4),'Median-', num2str(QQ),'-',W(1:length(W)-4),'.fig'];
% % saveas(hand,name);
% % % result(1,QQ)=psnrval;
% % % result(2,QQ)=ncval;
%  end


%  
% % % %  % % % % % % % 剪切
% % 
%  for QQ=2:2:4
%    wb=imread('watermarked.bmp');
%    wb(1:64*QQ,1:64*QQ,:)=0; % Cropping attack
%    imwrite(wb,'temp.bmp');
% 
%   for level=1:3
%     ws=double(wb(:,:,level));
%     coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
%      wa=coeffs{1,1};%%su3
%     wa=ws;
%     for  counter=1:lenw
%         i=position(1,counter);
%         j=position(2,counter);
% 
%         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%          h5 = @(r) block;
%      [s2,k2,r1,I2,K2,R2,h3]=dht(h5,R,N,n);
%      s2=s2*P;
%         Rext=fix(s2(1,1)/Qfac);
%         ext=mod(Rext,2);
%         if ext==0
%             watermark='0';
%         else
%             watermark='1';
%         end
%         ExWater(1,counter)=watermark;  
%      end
%      wateronedim(level,:)=ExWater;
% 
% end
%  total_bits = lenw; % 每个层次的总位数
% error_count = 0; % 错误位计数
% 
% for level = 1:3
%     for counter = 1:total_bits
%         original_bit = levelwatermark(level, counter);
%         extracted_bit = wateronedim(level, counter);
% 
%         if original_bit ~= extracted_bit
%             error_count = error_count + 1;
%         end
%     end
% end
% 
% ber = error_count / (total_bits * 3); % 总误码率
% disp(['误码率 (BER): ', num2str(ber)]);
% %水印回复
% for level=1:3
%     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
%     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% end
% 
% psnrval(k,m)=colorpsnr(wb,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(1+QQ);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
% subplot(122),imshow(uint8(wb)),title([num2str(QQ),'  PSNRr=',num2str(psnrval(k,m))]);
% name=[H(1:length(H)-4),'Crop_', num2str(QQ*0.125),'-',W(1:length(W)-4),'.fig'];
% saveas(hand,name);
% % % result(1,QQ)=psnrval;
% % % result(2,QQ)=ncval;
% m=m+1;
% end
% % 
% % 
% % % % watermarkedim=imread('watermarked.bmp');
% % % %  for QQ=1:5
% % % % 
% % % %   %h = fspecial('average');%('gaussian');%, 50, 45);%%('average')%,'disk','gaussian',
% % % % 
% % % %   for level=1:3
% % % %       %%wb = imfilter(watermarkedim(:,:,level), h);
% % % %       wb=medfilt2(watermarkedim(:,:,level),[QQ,1]);
% % % %        ttttt(:,:,level)=wb;
% % % %     ws=double(wb);
% % % %     waterposition=1;
% % % %      coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
% % % %      wa=coeffs{1,1};%%su3
% % % % for  counter=1:lenw
% % % %     i=position(1,counter);
% % % %     j=position(2,counter);
% % % % 
% % % %         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
% % % %          h5 = @(r) block;
% % % %      [s2,k2,r1,I2,K2,R2,h3]=dht(h5,R,N,n);
% % % %         Rext=fix(s2(1,1)/Qfac);
% % % %         ext=mod(Rext,2);
% % % %         if ext==0
% % % %             watermark='0';
% % % %         else
% % % %             watermark='1';
% % % %         end
% % % %         ExWater(1,counter)=watermark;  
% % % % end
% % % %      wateronedim(level,:)=ExWater;
% % % % 
% % % % end
% % % % 
% % % % %%%水印回复
% % % % 
% % % % for level=1:3
% % % %     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
% % % %     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% % % % 
% % % %  end
% % % % imwrite(ttttt,'temp.bmp');
% % % % psnrval=colorpsnr(ttttt,Host);
% % % % imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% % % % ncval=nc(uint8(extractlevelwatermark),uint8(Water));
% % % % hand=figure(QQ+2);
% % % % subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval)]);
% % % % subplot(122),imshow(uint8(ttttt)),title([num2str(QQ),'  PSNRr=',num2str(psnrval)]);
% % % % name=[H(1:length(H)-4),'Median-', num2str(QQ),'-',W(1:length(W)-4),'.fig'];
% % % % % saveas(hand,name);
% % % % % result(1,QQ)=psnrval;
% % % % % result(2,QQ)=ncval;
% % % %  end
% % 
% % 
% % % % % % % % % % 缩放攻击
%  a=imread('watermarked.bmp');
% for QQ=1
%    wb=imresize(a,QQ*0.5); 
%    wb=imresize(wb,[512 512]);
%    imwrite(wb,'temp.bmp');
%    for level=1:3
%     ws=double(wb(:,:,level));
%     waterposition=1;
%      % % coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
%      % % wa=coeffs{1,1};%%su3
%      wa=ws;
%  for  counter=1:lenw
%     i=position(1,counter);
%     j=position(2,counter);
%         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%        h5 = @(r) block;
%     [s2,k2,r1,I2,K2,R2,h3]=dht(h5,R,N,n);
%     s2=s2*P;
%         Rext=fix(s2(1,1)/Qfac);
%         ext=mod(Rext,2);
%         if ext==0
%             watermark='0';
%         else
%             watermark='1';
%         end
% 
%         ExWater(1,counter)=watermark;  
%   end
%        wateronedim(level,:)=ExWater;
% 
%  end
% % 
%  total_bits = lenw; % 每个层次的总位数
% error_count = 0; % 错误位计数
% 
% for level = 1:3
%     for counter = 1:total_bits
%         original_bit = levelwatermark(level, counter);
%         extracted_bit = wateronedim(level, counter);
% 
%         if original_bit ~= extracted_bit
%             error_count = error_count + 1;
%         end
%     end
% end
% 
% ber = error_count / (total_bits * 3); % 总误码率
% disp(['误码率 (BER): ', num2str(ber)]);
% %水印回复
% for level=1:3
%     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
%     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% end
% m=m+1;
% psnrval(k,m)=colorpsnr(wb,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(QQ+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
% subplot(122),imshow(uint8(wb)),title([num2str(QQ),'  PSNRr=',num2str(psnrval(k,m))]);
% name=[H(1:length(H)-4),'Resize_', num2str(QQ),'-',W(1:length(W)-4),'.fig'];psnrval(k,m)=colorpsnr(wb,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(QQ+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
% subplot(122),imshow(uint8(wb)),title([num2str(QQ),'  PSNRr=',num2str(psnrval(k,m))]);
% name=[H(1:length(H)-4),'Resize_', num2str(QQ*0.5),'-',W(1:length(W)-4),'.fig'];
% % saveas(hand,name);
% % result(1,QQ)=psnrval;
% % result(2,QQ)=ncval;
% end
% k=k+1;
% % % 
% % % 
% % % % % % 中值滤波
% % % % % 
% %   h = fspecial('average');%('gaussian');%, 50, 45);%%('average')%,'disk','gaussian',
% % 
% %  watermarkedim=imread('watermarked.bmp');
% %  for QQ=3:2:5
% % 
% %   h = fspecial('average');%('gaussian');%, 50, 45);%%('average')%,'disk','gaussian',
% % 
% %   for level=1:3
% %       %wb = imfilter(watermarkedim(:,:,level), h);
% %       wb=medfilt2(watermarkedim(:,:,level),[QQ,1]);
% %        ttttt(:,:,level)=wb;
% %     ws=double(wb);
% %     waterposition=1;
% %      % coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
% %      % wa=coeffs{1,1};%%su3
% %      wa = ws;
% % for  counter=1:lenw
% %     i=position(1,counter);
% %     j=position(2,counter);
% % 
% %         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
% %          h5 = @(r) block;
% %     [s1,k2,r1,I2,K2,R2,h3]=dht(h5,R,N,n);
% %     s1=s1*P;
% %         Rext=fix(s1(1,1)/Qfac);
% %         ext=mod(Rext,2);
% %         if ext==0
% %             watermark='0';
% %         else
% %             watermark='1';
% %         end
% %         ExWater(1,counter)=watermark;  
% % end
% %      wateronedim(level,:)=ExWater;
% % 
% % end
% % total_bits = lenw; % 每个层次的总位数
% % error_count = 0; % 错误位计数
% % 
% % for level = 1:3
% %     for counter = 1:total_bits
% %         original_bit = levelwatermark(level, counter);
% %         extracted_bit = wateronedim(level, counter);
% % 
% %         if original_bit ~= extracted_bit
% %             error_count = error_count + 1;
% %         end
% %     end
% % end
% % 
% % ber = error_count / (total_bits * 3); % 总误码率
% % disp(['误码率 (BER): ', num2str(ber)]);
% % %%水印回复
% % 
% % for level=1:3
% %     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
% %     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% % 
% % end
% %  m=m+1;
% % imwrite(ttttt,'temp.bmp');
% % psnrval(k,m)=colorpsnr(ttttt,Host);
% % imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% % ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
% % hand=figure(QQ+2);
% % subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
% % subplot(122),imshow(uint8(ttttt)),title([num2str(QQ),'  PSNRr=',num2str(psnrval(k,m))]);
% % name=[H(1:length(H)-4),'Median-', num2str(QQ),'-',W(1:length(W)-4),'.fig'];
% % saveas(hand,name);
% % % result(1,QQ)=psnrval;
% % % result(2,QQ)=ncval;
% %  end
% % % % 

% % 
% %旋转攻击
%  for Q=5
%                 watermarkedim=imread('watermarked.bmp');
%                 s1=imrotate(watermarkedim,Q);
%                 s2=imrotate(s1,-Q);
%                 s3=imresize(s2,[512 512]);
%                 attackedimg=RemovingBlackEdges(s3,512,512);
%                 figure(1);
%                 subplot(131);imshow(s1);
%                 subplot(132);imshow(s2);
%                 subplot(133);imshow(attackedimg);
%    for level=1:3
%     ws=double(attackedimg(:,:,level));
%                     waterposition=1;
%     waterposition=1;
%      coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
%      wa=coeffs{1,1};%%su3
%      wa=ws;
% for  counter=1:lenw
%     i=position(1,counter);
%     j=position(2,counter);
%         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%             h5 = @(r) block;
%     [s,k1,r,I1,K1,R1,h2]=dht(h5,R,N,n);
%      s=s*P;
%         Rext=fix(s(1,1)/Qfac);
%         ext=mod(Rext,2);
%         if ext==0
%             watermark='0';
%         else
%             watermark='1';
%         end
%         ExWater(1,counter)=watermark;  
%  end
%      wateronedim(level,:)=ExWater;
%  end
% total_bits = lenw; % 每个层次的总位数
% error_count = 0; % 错误位计数
% 
% for level = 1:3
%     for counter = 1:total_bits
%         original_bit = levelwatermark(level, counter);
%         extracted_bit = wateronedim(level, counter);
% 
%         if original_bit ~= extracted_bit
%             error_count = error_count + 1;
%         end
%     end
% end
% % 
% ber = error_count / (total_bits * 3); % 总误码率
% disp(['误码率 (BER): ', num2str(ber)]);
% % % % % % % 水印回复
% for level=1:3
%     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
%     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% 
% end
% 
% psnrval(k,1)=colorpsnr(attackedimg,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval(k,1)=nc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,1))]);
% subplot(122),imshow(uint8(attackedimg)),title([num2str(Q),'  PSNRr=',num2str(psnrval(k,1))]);
% name=[H(1:length(H)-4),'Rotaion5', num2str(Q),W(1:length(W)-4),'.fig'];
% saveas(hand,name);
% end
% % 

%%% % %%%%% Butterlow-pass filtering 
% 
% watermarkedim=imread('watermarked.bmp');
%  for QQ=6:4:10
%     for level=1:3
%         Ib=watermarkedim(:,:,level);
%         f=double(Ib);
%         g=fft2(f);
%         g=fftshift(g);
%         [N1,N2]=size(g);
%         d0=100;
%         n1=fix(N1/2);
%         n2=fix(N2/2);
%         for i=1:N1
%             for j=1:N2
%                 d=sqrt((i-n1)^2+(j-n2)^2);
%                 h=1/(1+0.414*(d/d0)^(2*QQ));
%                 result(i,j)=h*g(i,j);
%             end
%         end
%         result=ifftshift(result);
%         X2=ifft2(result);
%         X3=uint8(real(X2));
%         ttttt(:,:,level)=X3;
%         wb=double(X3);
%             % % % % % % % % Nonsubsampled Contourlet decomposition
%     % coeffs = nsctdec( wb, nlevels, dfilter, pfilter );
%     % wa=coeffs{1,1};%%su3
%     wa=wb;
%     for  counter=1:lenw
%     i=position(1,counter);
%     j=position(2,counter);
%         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%        h5 = @(r) block;
%     [s,k1,r,I1,K1,R1,h2]=dht(h5,R,N,n);
%     s=s*P;
%         Rext=fix(s(1,1)/Qfac);
%         ext=mod(Rext,2);
%         if ext==0
%             watermark='0';
%         else
%             watermark='1';
%         end
%         ExWater(1,counter)=watermark;  
%  end
%      wateronedim(level,:)=ExWater;
%     end
%     total_bits = lenw; % 每个层次的总位数
% error_count = 0; % 错误位计数
% 
% for level = 1:3
%     for counter = 1:total_bits
%         original_bit = levelwatermark(level, counter);
%         extracted_bit = wateronedim(level, counter);
% 
%         if original_bit ~= extracted_bit
%             error_count = error_count + 1;
%         end
%     end
% end
% 
% ber = error_count / (total_bits * 3); % 总误码率
% disp(['误码率 (BER): ', num2str(ber)]);
% % % 水印回复
% for level=1:3
%     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
%     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% 
% end
% imwrite(ttttt,'temp.bmp');
% psnrval(k,QQ)=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval(k,QQ)=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(QQ);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,QQ))]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(QQ),'  PSNRr=',num2str(psnrval(k,QQ))]);
% name=[H(1:length(H)-4),'ButterLowpass', '-',num2str(QQ),'-',W(1:length(W)-4),'.fig'];
% saveas(hand,name)
%  end
% k=k+1;
% 
% %%%%%%%%%%%%%%%%%guassianlowpass
% % % watermarkedim=imread('watermarked.bmp');
% % %  for QQ=6:6
% % % 
% % %   %h = fspecial('average');%('gaussian');%, 50, 45);%%('average')%,'disk','gaussian',
% % % filter=fspecial('gaussian',[QQ,QQ]);
% % %   for level=1:3
% % %       wb = imfilter(watermarkedim(:,:,level), filter);
% % % %       wb=medfilt2(watermarkedim(:,:,level),[QQ,1]);
% % %        ttttt(:,:,level)=wb;
% % %     ws=double(wb);
% % %     waterposition=1;
% % %      coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
% % %      wa=coeffs{1,1};%%su3
% % % for  counter=1:lenw
% % %     i=position(1,counter);
% % %     j=position(2,counter);
% % % 
% % %         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
% % %         h5 = @(r) block;
% % %    [s2,k1,r,I1,K1,R1,h2]=dht(h5,R,N,n);
% % %         Rext=fix(s2(1,1)/Qfac);
% % %         ext=mod(Rext,2);
% % %         if ext==0
% % %             watermark='0';
% % %         else
% % %             watermark='1';
% % %         end
% % %         ExWater(1,counter)=watermark;  
% % % end
% % %      wateronedim(level,:)=ExWater;
% % % 
% % % end
% % %    q=0;
% % %             for level=1:3
% % %                 for i=1:8192
% % %                     if wateronedim(level,i)~=levelwatermark(level,i)
% % %                         q=q+1;
% % %                     end
% % %                 end
% % %             end
% % %             ber=q/(8192*3);
% % % %%%水印回复
% % % for level=1:3
% % %     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
% % %     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% % % end
% % % 
% % % 
% % % imwrite(ttttt,'temp.bmp');
% % % psnrval=colorpsnr(ttttt,Host);
% % % imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% % % ncval=nc(uint8(extractlevelwatermark),uint8(Water));
% % % % % % % hand=figure(QQ+2);
% % % figure(4);
% % % subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc/ber=',num2str(ncval),'\',num2str(ber)]);
% % % subplot(122),imshow(uint8(ttttt)),title([num2str(QQ),'  PSNRr=',num2str(psnrval)]);
% % % % imwrite(uint8(extractlevelwatermark),'AAAGuass.jpg');
% % % % % name=[H(1:length(H)-4),'Median-', num2str(QQ),'-',W(1:length(W)-4),'.fig'];
% % % % saveas(hand,name);
% % % % % result(1,QQ)=psnrval;
% % % % % result(2,QQ)=ncval;
% % %  end
% % % 
% 
% % % %%平移攻击
% % % 
% % %                 watermarkedim=imread('watermarked.bmp');
% % % 
% % %                 attackedimg=TranslatingAttack(watermarkedim,20,40);
% % % %                 attackedimg=TranslatingAttack(attacked1,-35,25);
% % %                 figure(8);
% % % %                 subplot(131);imshow(s1);
% % % %                 subplot(132);imshow(s2);
% % % %                 subplot(133);imshow(attackedimg);
% % %    for level=1:3
% % %     ws=double(attackedimg(:,:,level));
% % %         % %             waterposition=1;
% % %     waterposition=1;
% % %      % coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
% % %      % wa=coeffs{1,1};%%su3
% % %      wa=ws;
% % % for  counter=1:lenw
% % %     i=position(1,counter);
% % %     j=position(2,counter);
% % %         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
% % %           h5 = @(r) block;
% % %      [s2,k1,r,I1,K1,R1,h2]=dht(h5,R,N,n);
% % %      s2=s2*P;
% % %         Rext=fix(s2(1,1)/Qfac);
% % %         ext=mod(Rext,2);
% % %         if ext==0
% % %             watermark='0';
% % %         else
% % %             watermark='1';
% % %         end
% % %         ExWater(1,counter)=watermark;  
% % %  end
% % %      wateronedim(level,:)=ExWater;
% % %  end
% % % q=0;
% % %             for level=1:3
% % %                 for i=1:8192
% % %                     if wateronedim(level,i)~=levelwatermark(level,i)
% % %                         q=q+1;
% % %                     end
% % %                 end
% % %             end
% % %             ber=q/(8192*3);
% % % 
% % % 
% % % %%水印回复
% % % for level=1:3
% % %     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
% % %     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来    
% % % end
% % % 
% % % psnrval(k,1)=colorpsnr(attackedimg,Host);
% % % imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% % % ncval(k,1)=nc(uint8(extractlevelwatermark),uint8(Water));
% % % % hand=figure(Q+2);
% % % subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc/ber=',num2str(ncval(k,1)),'/',num2str(ber)]);
% % % subplot(122),imshow(uint8(attackedimg));
% % % % name=[H(1:length(H)-4),'Rotaion5', num2str(Q),W(1:length(W)-4),'.fig'];
% % % % saveas(hand,name);
% % % % result(1,QQ)=psnrval;
% % % % result(2,QQ)=ncval;
% % % imwrite(uint8(extractlevelwatermark),'AAAtranslating.jpg');
% % % 
% 
% 
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Affine仿射
% %  for para=0.09
% %      watermarkedim=imread('watermarked.bmp');
% %  attackedimg=double(AffiningAttack(watermarkedim,1,para,para)); 
% % %                         figure(1);
% % %                 subplot(131);imshow(s1);
% % %                 subplot(132);imshow(s2);
% % %                 subplot(133);imshow(attackedimg);
% % % imshow(uint8(attacked));
% %    for level=1:3
% %     ws=double(attackedimg(:,:,level));
% %         % %             waterposition=1;
% %     waterposition=1;
% %      % coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
% %      % wa=coeffs{1,1};%%su3
% %      wa = ws;
% % for  counter=1:lenw
% %     i=position(1,counter);
% %     j=position(2,counter);
% %         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
% %        h5 = @(r) block;
% %      [s2,k1,r,I1,K1,R1,h2]=dht(h5,R,N,n);
% %         Rext=fix(s2(1,1)/Qfac);
% %         ext=mod(Rext,2);
% %         if ext==0
% %             watermark='0';
% %         else
% %             watermark='1';
% %         end
% %         ExWater(1,counter)=watermark;  
% %  end
% %      wateronedim(level,:)=ExWater;
% %  end
% %  q=0;
% %             for level=1:3
% %                 for i=1:8192
% %                     if wateronedim(level,i)~=levelwatermark(level,i)
% %                         q=q+1;
% %                     end
% %                 end
% %             end
% %             ber=q/(8192*3);
% % %%水印回复
% % for level=1:3
% %     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
% %     extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
% % 
% % end
% % 
% % % psnrval(k,1)=colorpsnr(attackedimg,Host);
% % imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% % ncval(k,1)=nc(uint8(extractlevelwatermark),uint8(Water));
% % % % hand=figure(Q+2);
% % figure(7);
% % subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc/ber=',num2str(ncval(k,1)),'/',num2str(ber)]);
% % subplot(122),imshow(uint8(attackedimg)),title([num2str(para),'  PSNRr=',num2str(psnrval(k,1))]);
% % 
% % end