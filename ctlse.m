clc;
clear all;
blocksize=4;
sgm=30;
% %按照由左至右，自上而下的顺序来间隔选择嵌入块 并采用类似DC的量化方法
H='Barbara.jpg';
W='ldu3232.jpg';
%W='200100732.jpg';
Host=imread(H);%读入宿主图片
Water=imread(W);%读入水印
lenw=size(Water,1)*size(Water,2)*8;%水印矩阵行数*列数，32*32*8水印图像所占的bit位吗
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
% Hostlevelback=Hostlevel;
nlevels =1;%[0, 1, 3] ;        % Decomposition level
pfilter = 'maxflat' ;%'9-7';%              % Pyramidal filter
dfilter = 'dmaxflat7' ; %  'pkva';           % Directional filter

% Nonsubsampled Contourlet decomposition
coeffs = nsctdec( Hostlevel, nlevels, dfilter, pfilter );
Lowpass=coeffs{1,1};%%save in su1

Waterold=levelwatermark(level,:);%将对应的水印level二进制序列赋给waterold，水印分了三层，一层一层的嵌入
[Hr,Hc]=size(Lowpass);%返回Hostlevel二维矩阵的行列数。
waterposition=1;%嵌入水印位置，初始值为1
for i=1:Hr/blocksize%将载体图像分块，按行可以分成几个。
        if mod(i,2)==1 && mod(level,2)==1%根据层数和行数构造一个网格选块布局
            star=1;
        else
            star=2;
        end
    for j=star:2:Hc/blocksize%先对列进行划分，列按照4*4分块原则可以分多少次，间隔选列，交叉进行，步长为2
        watermark=Waterold(1,waterposition);%将此时该嵌入的水印的数值赋给水印位
        block=Lowpass((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);%分块位置
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
        if (watermark=='1') %水印嵌入值等于1的时候，量化过程如下。
           if mod(lmt+1,2)==1
              s1=(lmt+1.5)*sgm1;
           else
                s1=(lmt+0.5)*sgm1;
           end
        else %水印嵌入值等于0的时侯
           if mod(lmt,2)==1
              s1=(lmt+1.5)*sgm1;
           else
               s1=(lmt+0.5)*sgm1;
           end
        end
        s(1,1)=s1;
        blocknew=u*s*u';
        Lowpass((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize)=blocknew;%将分块值重新赋予原来的分块
        waterposition=waterposition+1;%水印嵌入位置（嵌入量）加1
        if waterposition>lenw %%%==2799%%判断是否将所有的水印序列都嵌入进去
            break;
        end
    end
     if waterposition>lenw %%==5778==5578==2799
           break;
     end
    end%这一层水印图像嵌入完成
    coeffs{1,1}=Lowpass;
    Hostlevel= nsctrec(coeffs, dfilter, pfilter) ;
watermarkedim(:,:,level)=Hostlevel;%将这一层的载体图像二维矩阵赋给
imwrite(uint8(Hostlevel),'www.bmp');%将Hostlevel（嵌入水印后的宿主图像）（原来为double类型的）转换成uint8,保存图片.png格式
imwrite(Host(:,:,level),'hhh.bmp');%保存宿主这一层的图像
psnrval=psnrval+PSNR('www.bmp','hhh.bmp');%将这三层加一块，计算两者图片之间的相似度
 end%RGB三层图像均已嵌入完成
hfx=watermarkedim;
imwrite(uint8(watermarkedim),'watermarked.bmp');%保存嵌入水印的图片为watermarked.bmp（后期攻击所用图片）
%psnrval=PSNR(H,'watermarked.bmp');
hand=figure(2),subplot(131),imshow(Host),title('Original image');
subplot(132),imshow('watermarked.bmp'),title([' unattack PSNR=',num2str(psnrval/3)]);%显示嵌入完成后的水印图片
toc;
%--------------------------------------------------
%水印提取
for level=1:3
    wa=watermarkedim(:,:,level);
    Waterold=levelwatermark(level,:);
    coeffs = nsctdec( wa, nlevels, dfilter, pfilter );
    wa=coeffs{1,1};%%su3
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

%%算法测试
%%
%  watermarkedim=imread('watermarked.bmp');
%  for Q=10:10:100
%     for level=1:3
%     imwrite(watermarkedim(:,:,level),'temp.bmp','jpg','quality',Q);%进行JPEG压缩攻击
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
%         
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
%  end
%  
% %%水印恢复
% 
% for level=1:3
%    extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
%     extractlevelwatermark(:,:,level)=Arnold(extractlevelwatermark(:,:,level),6,1);
%  end
% imwrite(ttttt,'temp.bmp');
% psnrval=colorpsnr(ttttt,Host);
% imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
% ncval=colornc(uint8(extractlevelwatermark),uint8(Water));
% hand=figure(Q/10+2);
% subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark JPEG nc=',num2str(ncval)]);
% subplot(122),imshow(uint8(ttttt)),title([num2str(Q),'  PSNRr=',num2str(psnrval)]);
% name=['lenaJPEG', num2str(Q/10),'307','.fig'];
% % saveas(hand,name)
% % result(ss,3+Q/10)=ncval;
% end
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
% %水印回复
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
%%噪声攻击
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
end
