 clc
clear all;
blocksize=4;
load position;
im={'house.jpg'};
k=1; 
for Qfac=26    %2:2:34
  m=1;
%Qfac=4;
% T=0.012 
%实验结果表明：block=4×4具有很强的不可见性，
%优越性明显；block=8×8具有较好的鲁棒性，效果与前者相差无几，
%在水印容量较大时可以考虑用较小的分块尺寸得到多个分块用来隐藏水印。 
%改变分别改变U和
%H='lena.jpg';%;'peppers.jpg';%'baboon.jpg';%'F16.jpg';%
%H=im{1,imm};
H='house.jpg';
W='3073232.tif';
Host=imread(H);
Water=imread(W);
lenw=size(Water,1)*size(Water,2)*8;
figure(1),subplot(1,2,1),imshow(Host),title('Host image');
subplot(1,2,2),imshow(Water),title('watermark image');

 t1=cputime;
  for level=1:3
     temp1=Water(:,:,level);%矩阵分层提取
%      figure(2),subplot(1,3,level),imshow(temp1);
     AfterHundun= Hundun(temp1,1);
     levelwatermark(level,:)=gainlevelwatermark(AfterHundun);
%      figure(3),subplot(1,3,level),imshow(AfterHundun);
  end
 


%水印嵌入
% t1=0;t2=0;t1=cputime;
for level=1:3
%     Hostimage=Host(:,:,level); %主图像降维
%     figure(4),subplot(1,3,level),imshow(Hostimage);
    Hostlevel=double(Host(:,:,level));
    Hostlevelback=Hostlevel;
    nlevels =1;%[0, 1, 3] ;        % Decomposition level
    pfilter = 'maxflat' ;%'9-7';%              % Pyramidal filter
    dfilter = 'dmaxflat7' ; %  'pkva';           % Directional filter

% Nonsubsampled Contourlet decomposition
coeffs = nsctdec( Hostlevel, nlevels, dfilter, pfilter );
Hostlevel=coeffs{1,1};%%save in su1

Waterold=levelwatermark(level,:);
[Hr,Hc]=size(Hostlevel);
waterposition=1;

for counter=1:lenw
    i=position(1,counter);%根据position取出应该的嵌入位置
    j=position(2,counter);
    watermark=Waterold(1,counter);%%根据水印序列取出水印位
    block=Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
    %将经过变换后的
    [u,s]=schur(block);
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
        blocknew=u*s*u';
    Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize)=blocknew;  
 end
% Reconstruct image
    coeffs{1,1}=Hostlevel;%%su2
    Hostlevel= nsctrec(coeffs, dfilter, pfilter) ;
    watermarkedim(:,:,level)=Hostlevel;
    imwrite(uint8(Hostlevel),'www.bmp');
    imwrite(Host(:,:,level),'hhh.bmp');
    ssim(level)=colorssim(watermarkedim(:,:,level),Host(:,:,level));
end
%  figure(6); subplot(121),imshow('www.bmp') ;subplot(122),imshow('hhh.bmp') ;
%result(k,1)=cputime-t1;
imwrite(uint8(watermarkedim),'watermarked.bmp');
psnrval(k,1)=colorpsnr(watermarkedim,Host);
ssimval(k,1)=(ssim(1)+ssim(2)+ssim(3))/3;
hand=figure(5);
subplot(131),imshow(Host),title('Original image');
subplot(132),imshow('watermarked.bmp'),title([' PSNR=',num2str(psnrval(k,1)),' SSIM=',num2str(ssimval(k,1)),' Qfac=',num2str(Qfac)]);

t2=cputime;
t3=t2-t1;
fprintf('The embedding time is %5f\n',t3);


 t4=cputime;
%水印提取
watermarkedim=double(imread('watermarked.bmp'));
for level=1:3
    wa=watermarkedim(:,:,level);
    Waterold=levelwatermark(level,:);
% Nonsubsampled Contourlet decomposition
    coeffs = nsctdec( wa, nlevels, dfilter, pfilter );
    wa=coeffs{1,1};%%su3
    for counter=1:lenw
        i=position(1,counter);
        j=position(2,counter);
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u1,s1]=schur(block);
        Rext=fix(s1(1,1)/Qfac);
        ext=mod(Rext,2);
        if ext==0
            watermark='0';
        else
            watermark='1';
        end
        ExWater(1,counter)=watermark;     
    end
    wateronedim(level,:)=ExWater;
end
%  result(k,2)=cputime-t2;
%--------------------------------------------------
%水印恢复
for level=1:3
    extractlevelwatermark (:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
    extractlevelwatermark(:,:,level)=Hundun(extractlevelwatermark(:,:,level),2);%将置乱后的像素值利用Arnold变换转变回来
 end
%extractlevelwatermark=watermarkrestore(wateronedim);
%toc;
imwrite(uint8(extractlevelwatermark),'extrwater.bmp');
ncval(k,1)=nc(uint8(extractlevelwatermark),uint8(Water));
subplot(133),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,1))]);

t5=cputime;
t6=t5-t4;
fprintf('The Extraciton time is %5f',t6);

%-----------------------------------------------------------------------------------------------------------------------------
% % JPEG2000 测试 一  图像整体攻击提取效果较差
%  watermarkedim=imread否认 难比 vc差不多你付出   从踩踩踩
% ('watermarked.bmp');
%  for QA=1:3:7
%   for level=1:3
%       imwrite(watermarkedim(:,:,level),'temp.bmp','jp2','compressionratio',QA);
%       wb=double(imread('temp.bmp'));
%       ttttt(:,:,level)=wb;
%     %Nonsubsampled Contourlet decomposition
%     coeffs = nsctdec( wb, nlevels, dfilter, pfilter );
%      wa=coeffs{1,1};%%su3
%     for  counter=1:lenw
%         i=position(1,counter);
%         j=position(2,counter);
%         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%          [u2,s2]=schur(block);
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
% %水印回复
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
% %xlswrite('D:\test\text.xls',num2str(ncval),'tt','C1:C4')
% % saveas(hand,name);
% % result(1,QA)=psnrval;
% % result(2,QA)=ncval;
% 
% end





end