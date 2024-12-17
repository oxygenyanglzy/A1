clc;
clear all;
blocksize=4;

% %�����������ң����϶��µ�˳�������ѡ��Ƕ��� ����������DC����������
ss=0;
for delta=38:38
   ss=ss+1;
   result(ss,1)=delta;%��delta��ֵ�����������ֻ��߽о���result�ĵ�һ��λ��
H='lena.jpg';
W='3073232.tif';
%W='200100732.jpg';
Host=imread(H);%��������ͼƬ
Water=imread(W);%����ˮӡ
lenw=size(Water,1)*size(Water,2)*8;%ˮӡ��������*������32*32*8ÿˮӡͼ����ռ��bitλ
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
Hostlevelback=Hostlevel;
Waterold=levelwatermark(level,:);%����Ӧ��ˮӡlevel���������и���waterold��ˮӡ�������㣬һ��һ���Ƕ��
[Hr,Hc]=size(Hostlevel);%����Hostlevel��ά�������������
waterposition=1;%Ƕ��ˮӡλ�ã���ʼֵΪ1
for i=1:Hr/blocksize%������ͼ��ֿ飬���п��Էֳɼ�����
        if mod(i,2)==1 && mod(level,2)==1%�⼸��������Ϊ���ָ����ѡ�飿��ͬʱ���������У�������ʱ��
            star=1;
        else
            star=2;
        end
    for j=star:2:Hc/blocksize%�ȶ��н��л��֣��а���4*4�ֿ�ԭ����Էֶ��ٴΣ����ѡ�У�������У�����Ϊ2
        watermark=Waterold(1,waterposition);%����ʱ��Ƕ���ˮӡ�ĸ�������ˮӡλ
        block=Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);%�ֿ�λ��ѡ��
        [u,r]=qr(double(block));%����QR�ֽ⣬
        %������QR�ֽ���U�ĵ�һ�ж��ǷǸ��������ԾͲ�������ֵ���⣬��SVD��SChur�ĵ�һ���Ǹ����Ϳ��Ǿ���ֵ����
        c=r(1,3);%����r(1,3)λ��Ƕ�����
        k=floor(ceil(c/delta)/2);%ceil������Ҫ������ȡ���������͸��������ԣ�y = floor(x) ������x��Ԫ��ȡ����ֵyΪ�����ڱ����������������ڸ������ֱ��ʵ�����鲿ȡ��
        if (watermark=='1') %ˮӡǶ��ֵ����1��ʱ�������������¡�
            c11=2*k*delta+0.5*delta;
            c12=2*(k-1)*delta+0.5*delta;
            cnew=c11;
            if abs(c-c12)<abs(c-c11)
                cnew=c12;
            end
        else %ˮӡǶ��ֵ����0��ʱ��
            c01=2*k*delta-0.5*delta;
            c02=2*(k+1)*delta-0.5*delta;
            cnew=c01;
            if abs(c-c02)<abs(c-c01)
                 cnew=c02;
            end
        end
        r(1,3)=cnew;%R����1��3λ�����¸�ֵ��
        blocknew=u*r;%��qr�ֽ�
        Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize)=blocknew;%���ֿ�ֵ���¸���ԭ���ķֿ�
         waterposition=waterposition+1%ˮӡǶ��λ�ã�Ƕ��������1
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
sqt=watermarkedim;
imwrite(uint8(watermarkedim),'watermarked.bmp');%����Ƕ��ˮӡ��ͼƬΪwatermarked.bmp�����ڹ�������ͼƬ��
%psnrval=PSNR(H,'watermarked.bmp');
hand=figure(2),subplot(131),imshow(Host),title('Original image');
subplot(132),imshow('watermarked.bmp'),title([' unattack PSNR=',num2str(psnrval/3)]);%��ʾǶ����ɺ��ˮӡͼƬ
toc;
%--------------------------------------------------
%ˮӡ��ȡ
watermarkedim=double(imread('watermarked.bmp'));
watermarkedim=sqt;%��ˮӡ��ͼ��
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
        tt=mod(ceil(ctnew/delta),2);%������ż�����ж�ͼ��Ƕ�����0����1
        if  tt==0   %
           watermark='0';
        else
            watermark='1';
        end
        
        
        ExWater(1,waterposition)=watermark;%����ȡ�����Ķ�����ֵ������ȡ������ ExWater����
        waterposition=waterposition+1;
        if waterposition>lenw %==5578%==2799
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
saveas(hand,'lena unattack results');
toc;
result(ss,2)=psnrval/3;
result(ss,3)=ncval;

% %--------------------------------------------------
% %�㷨����
% �����Ƕ�ͼ����и��ֹ������ԣ����Ƿ�����ȡ��ˮӡ
% %---------------------------------------------------------------------
%JPEG ���� ��  ͼ��ֲ㹥����ȡЧ���Ϻ�

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
 
%ˮӡ�ظ�

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
% % % %��������  Ч��Ҫ��һЩ��������4��
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
 
% ˮӡ�ظ�

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
 
%ˮӡ�ظ�

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
% % Median filtering Ч��Ҫ��һЩ��������4��
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
 
%ˮӡ�ظ�

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
% % %%% Butterlow-pass filtering Ч��Ҫ��һЩ��������4��
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
 %ˮӡ�ظ�
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

% % %resize����  Ч��Ҫ�����Ϻã�������4��
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
 
%ˮӡ�ظ�

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
 % %resize����  Ч��Ҫ�����Ϻã�������4��
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
 
%ˮӡ�ظ�

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
% % % ���й��� 
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
 
% % ˮӡ�ظ�

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
% % % %rotation����  ��ת������Ч�����ã�
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
% %ˮӡ�ظ�
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
% %    %%%%-------blurring����  Ч��Ҫ������������4��
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
