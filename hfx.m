clc;
clear all;
blocksize=4;

% %按照由左至右，自上而下的顺序来间隔选择嵌入块 并采用类似DC的量化方法
ss=0;
for delta=38:38
   ss=ss+1;
   result(ss,1)=delta;%将delta的值赋给变量（又或者叫矩阵）result的第一个位置
H='lena.jpg';
W='3073232.tif';
%W='200100732.jpg';
Host=imread(H);%读入宿主图片
Water=imread(W);%读入水印
lenw=size(Water,1)*size(Water,2)*8;%水印矩阵行数*列数，32*32*8每水印图像所占的bit位
figure(1),subplot(121),imshow(Host),title('Host image');
subplot(122),imshow(Water),title('watermark image');%显示水印图片和宿主图片
%--------------------------------------------------
%原始水印处理 调用子函数形成分层水印levelwatermark
for level=1:3
    temp1=Water(:,:,level);%矩阵分层提取
    afterarnold=Arnold(temp1,6,0);%将水印图像分层置乱，并迭代6次，加密
    levelwatermark(level,:)=gainlevelwatermark(afterarnold);%获得加密后的分层图像矩阵，并转换成二进制数据流水印信息
end
%--------------------------------------------------
%水印嵌入
tic;
psnrval=0;
 for level=1:3
Hostlevel=double(Host(:,:,level));%先将（载体图像）分层并转换成double类型的二维矩阵。
Hostlevelback=Hostlevel;
Waterold=levelwatermark(level,:);%将对应的水印level二进制序列赋给waterold，水印分了三层，一层一层的嵌入
[Hr,Hc]=size(Hostlevel);%返回Hostlevel二维矩阵的行列数。
waterposition=1;%嵌入水印位置，初始值为1
for i=1:Hr/blocksize%将载体图像分块，按行可以分成几个。
        if mod(i,2)==1 && mod(level,2)==1%这几行意欲何为？分隔间隔选块？（同时满足奇数行，奇数层时）
            star=1;
        else
            star=2;
        end
    for j=star:2:Hc/blocksize%先对列进行划分，列按照4*4分块原则可以分多少次，间隔选列，交叉进行，步长为2
        watermark=Waterold(1,waterposition);%将此时该嵌入的水印的个数赋给水印位
        block=Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);%分块位置选择
        [u,r]=qr(double(block));%进行QR分解，
        %由于在QR分解中U的第一列都是非负数，所以就不考绝对值问题，而SVD或SChur的第一列是负数就考虑绝对值问题
        c=r(1,3);%测试r(1,3)位置嵌入情况
        k=floor(ceil(c/delta)/2);%ceil函数主要是向上取整，整数和复数都可以，y = floor(x) 函数将x中元素取整，值y为不大于本身的最大整数。对于复数，分别对实部和虚部取整
        if (watermark=='1') %水印嵌入值等于1的时候，量化过程如下。
            c11=2*k*delta+0.5*delta;
            c12=2*(k-1)*delta+0.5*delta;
            cnew=c11;
            if abs(c-c12)<abs(c-c11)
                cnew=c12;
            end
        else %水印嵌入值等于0的时候
            c01=2*k*delta-0.5*delta;
            c02=2*(k+1)*delta-0.5*delta;
            cnew=c01;
            if abs(c-c02)<abs(c-c01)
                 cnew=c02;
            end
        end
        r(1,3)=cnew;%R矩阵1，3位置重新赋值。
        blocknew=u*r;%逆qr分解
        Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize)=blocknew;%将分块值重新赋予原来的分块
         waterposition=waterposition+1%水印嵌入位置（嵌入量）加1
        if waterposition>lenw %%%==2799%%判断是否将所有的水印序列都嵌入进去
            break;
        end
    end
     if waterposition>lenw %%==5778==5578==2799
           break;
     end
    end%这一层水印图像嵌入完成
watermarkedim(:,:,level)=Hostlevel;%将这一层的载体图像二维矩阵赋给
imwrite(uint8(Hostlevel),'www.bmp');%将Hostlevel（嵌入水印后的宿主图像）（原来为double类型的）转换成uint8,保存图片.png格式
imwrite(Host(:,:,level),'hhh.bmp');%保存宿主这一层的图像
psnrval=psnrval+PSNR('www.bmp','hhh.bmp');%将这三层加一块，计算两者图片之间的相似度
 end%RGB三层图像均已嵌入完成
sqt=watermarkedim;
imwrite(uint8(watermarkedim),'watermarked.bmp');%保存嵌入水印的图片为watermarked.bmp（后期攻击所用图片）
%psnrval=PSNR(H,'watermarked.bmp');
hand=figure(2),subplot(131),imshow(Host),title('Original image');
subplot(132),imshow('watermarked.bmp'),title([' unattack PSNR=',num2str(psnrval/3)]);%显示嵌入完成后的水印图片
toc;
%--------------------------------------------------
%水印提取
watermarkedim=double(imread('watermarked.bmp'));
watermarkedim=sqt;%带水印的图像
tic;
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
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u,r]=qr(double(block));
        ctnew=r(1,3);
        tt=mod(ceil(ctnew/delta),2);%根据奇偶性来判断图像嵌入的是0还是1
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
        end
        
        
        ExWater(1,waterposition)=watermark;%将提取出来的二进制值赋给提取出来的 ExWater矩阵
        waterposition=waterposition+1;
        if waterposition>lenw %==5578%==2799
            break;
        end
              
    end
     if waterposition>lenw %%==5578
           break;
     end
    end
     wateronedim(level,:)=ExWater;%将相应的水印值赋给相应的R,G,B层，
    end

% %--------------------------------------------------
%水印恢复

for level=1:3
    extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));%将水印二进制序列转换成相应的像数值，
    extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);%将置乱后的像素值利用Arnold变换转变回来
    
 end
% wateronedim=reshape(ExWater,size(Water,1),size(Water,2));

imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=nc(uint8(extractlevelwatermark),uint8(Water));%比较水印图像nc值
subplot(133),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval)]);
saveas(hand,'lena unattack results');
toc;
result(ss,2)=psnrval/3;
result(ss,3)=ncval;

% %--------------------------------------------------
% %算法测试
% 以下是对图像进行各种攻击测试，看是否能提取出水印
% %---------------------------------------------------------------------
%JPEG 测试 二  图像分层攻击提取效果较好

 watermarkedim=imread('watermarked.bmp');
 for Q=10:10:100
    for level=1:3
    imwrite(watermarkedim(:,:,level),'temp.bmp','jpg','quality',Q);%进行JPEG压缩攻击
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
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
       
         [u,r]=qr(double(block));
        ctnew=r(1,3);
        tt=mod(ceil(ctnew/delta),2);
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
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
 
%%水印恢复

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
saveas(hand,name)
result(ss,3+Q/10)=ncval;
end
 ss2=3+Q/10;
% % 
% %JPEG2000
 watermarkedim=imread('watermarked.bmp');
 for Q=1:10
    for level=1:3
        imwrite(watermarkedim(:,:,level),'temp.bmp','jp2','compressionratio',Q);
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
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
         [u,r]=qr(double(block));
         ctnew=r(1,3);
        tt=mod(ceil(ctnew/delta),2);
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
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
 
%水印回复

for level=1:3
     extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
    extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
hand=figure(Q+2);
subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark JPEG2000 nc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenaJPEG2000', num2str(Q),'307','.fig'];
saveas(hand,name)

result(ss,ss2+Q)=ncval;
 end
ss2=ss2+Q;
% 
% 
% % 
% % % %---------------------------------------------------------------------
% % % 
% % % %噪声攻击  效果要好一些（比文献4）
watermarkedim=imread('watermarked.bmp');
 for Q=1:5
    for level=1:3
             wb=imnoise(watermarkedim(:,:,level),'gaussian',0,0.001*Q);
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
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u,r]=qr(double(block));
        ctnew=r(1,3);
        tt=mod(ceil(ctnew/delta),2);
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
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
 
% 水印回复

for level=1:3
    extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
    extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
hand=figure(Q+2);
subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark gaussiannc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenagaussian', num2str(Q),'by 0.1 right','.fig'];
saveas(hand,name)
result(ss,ss2+Q)=ncval;
end
ss2=ss2+Q;

%  % % -salt noise----
 watermarkedim=imread('watermarked.bmp');
 for Q=1:5
    for level=1:3
       wb=imnoise(watermarkedim(:,:,level),'salt & pepper',0.02*Q);  %salt & pepper
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
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
         [uu,rr]=qr(block);
         ctnew=rr(1,3);
        tt=mod(ceil(ctnew/delta),2);
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
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
 
%水印回复

for level=1:3
    extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
    extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
hand=figure(Q+2);
subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark salt nc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenaSalt', num2str(Q),'by0.02','.fig'];
saveas(hand,name)
result(ss,ss2+Q)=ncval;
 end

% % %---------------------------------------------------------------------
% % 
% % Median filtering 效果要好一些（比文献4）
% 
ss2=ss2+Q;
watermarkedim=imread('watermarked.bmp');
 for Q=1:5
    for level=1:3
      wb=medfilt2(watermarkedim(:,:,level),[Q,Q]);
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
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u,r]=qr(block);

        ctnew=r(1,3);
        tt=mod(ceil(ctnew/delta),2);
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
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
 
%水印回复

for level=1:3
      extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
    extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);;
end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
hand=figure(Q+2);
subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark Median nc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenaMedian', num2str(Q),'by','.fig'];
saveas(hand,name)
result(ss,ss2+Q)=ncval;
 end
ss2=ss2+Q;
% % %%% Butterlow-pass filtering 效果要好一些（比文献4）
% 
watermarkedim=imread('watermarked.bmp');
 for Q=1:5
    for level=1:3
   Ib=watermarkedim(:,:,level);
f=double(Ib);
g=fft2(f);
g=fftshift(g);
[N1,N2]=size(g);
d0=100;
n1=fix(N1/2);
n2=fix(N2/2);
for i=1:N1
    for j=1:N2
        d=sqrt((i-n1)^2+(j-n2)^2);
        h=1/(1+0.414*(d/d0)^(2*Q));
        resul(i,j)=h*g(i,j);
         end
end
resul=ifftshift(resul);
X2=ifft2(resul);
X3=uint8(real(X2));
    ttttt(:,:,level)=X3;
    wa=double(X3);
    waterposition=1;
    for i=1:Hr/blocksize
       if mod(i,2)==1 && mod(level,2)==1
            star=1;
        else
            star=2;
        end
    for j=star:2:Hc/blocksize
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u,r]=qr(block);
         ctnew=r(1,3);
        tt=mod(ceil(ctnew/delta),2);
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
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
 %水印回复
for level=1:3
 extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
   extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
hand=figure(Q+2);
subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark Butterlow-pass nc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenaButterLowpass', 'by',num2str(Q),'.fig'];
saveas(hand,name)
result(ss,ss2+Q)=ncval;
end

% % %resize攻击  效果要调整较好（比文献4）
ss2=ss2+Q;
watermarkedim=imread('watermarked.bmp');
 for Q=2:5
    for level=1:3
        wb=imresize(watermarkedim(:,:,level),1/Q); 
        wb=imresize(wb,[512 512]);   
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
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u,r]=qr(block);
         ctnew=r(1,3);
        tt=mod(ceil(ctnew/delta),2);
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
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
 
%水印回复

for level=1:3
 extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
hand=figure(Q+2);
subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark Resize small nc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenaSCaling','by 111', num2str(Q),'.fig'];
saveas(hand,name)
result(ss,ss2+Q)=ncval;
 end

 ss2=ss2+Q;
 % %resize攻击  效果要调整较好（比文献4）
watermarkedim=imread('watermarked.bmp');
 for Q=2:5
    for level=1:3
        wb=imresize(watermarkedim(:,:,level),Q); %gausian
        wb=imresize(wb,[512 512]);   
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
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u,r]=qr(block);
         ctnew=r(1,3);
        tt=mod(ceil(ctnew/delta),2);
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
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
 
%水印回复

for level=1:3
 extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
hand=figure(Q+2);
subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark ResizeBig nc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenaSCaling','times ', num2str(Q),'.fig'];
saveas(hand,name)
result(ss,ss2+Q)=ncval;
 end

ss2=ss2+Q;
% % %     %wb(100:200,:,:)=255;% Cropping attack
% % % 剪切攻击 
watermarkedim=imread('watermarked.bmp');
wb=watermarkedim;
 for Q=1:4
   wb(1:64*Q,:,:)=255;
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
        block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        [u,r]=qr(block);
         ctnew=r(1,3);
        tt=mod(ceil(ctnew/delta),2);
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
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
 
% % 水印回复

for level=1:3
    tt=gainlevelimage(wateronedim(level,:));
    te=Arnold(tt,6,1);
    extimage(:,:,level)=te;
    imshow(te);
    
 end
imwrite(ttttt,'temp.bmp');
psnrval=colorpsnr(ttttt,Host);
imwrite(uint8( extimage),'extrwater.bmp')
ncval=colornc(uint8( extimage),uint8(Water));
hand=figure(Q+2);
subplot(121),imshow( uint8( extimage)),title(['Extracted Watermark Crop nc=',num2str(ncval)]);
subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
name=['lenaCrop64by',num2str(Q), '.fig'];
saveas(hand,name)
result(ss,ss2+Q)=ncval;
 end
 end

% 
% 
% 
% % % %rotation攻击  旋转攻击的效果不好）
% a=imread('lenaRotation30end.bmp');
% 
%  for Q=1:1
%      
%      wb=a; %gausian
% %    wb=imrotate(wb,-Q);
% %    wb=imrotate(wb,-Q*0.1); %gausian
%    imshow(wb);
%    wb=imresize(wb,[512 512]);
%     for level=1:3
%     wa=double(wb(:,:,level));
%      ttttt=double(wb(:,:,level));
%     waterposition=1;
%     for i=1:Hr/blocksize
%        for j=1:Hc/blocksize
%         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%         [u,s]=qr(block);
%         if abs( u(2,1))< abs(u(3,1))
%            watermark='0';
%         else
%             watermark='1';
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
% %水印回复
% 
% for level=1:3
%  extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
%     extractlevelwatermark(:,:,level)=arnold(extractlevelwatermark(:,:,level),6,1);
% end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(wb,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenaRotation', num2str(Q),'by 30','.fig'];
% saveas(hand,name)
% result(1,Q)=psnrval;
% result(2,Q)=ncval;
%  end
%  
% 
%  
% %    %%%%-------blurring攻击  效果要调整（比文献4）
% %    wb=imread('lenarotation5end.bmp');%blurring01
% %     %wb=imread('lenablur02.bmp');%blurring1
% %        %wb=imread('lenasharp01.bmp');%Sharpen0.1
% %     % wb=imread('lenasharp02.bmp');%Sharpen1.0
% %    %wb=imread('lenacontrast.bmp');%%contrastlenadefault
% %    
% % % 
%   watermarkedim=imread('watermarked.bmp');
%   wb=watermarkedim;  
%   wb(256:257,256:257,:)=255;

% 
% 
% 
