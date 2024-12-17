clc;
clear all;
blocksize=4;

% %按照由左至右，自上而下的顺序来间隔选择嵌入块 并采用类似DC的量化方法
ss=0;
% for delta=38:38
%    ss=ss+1;
%    result(ss,1)=delta;%将delta的值赋给变量（又或者叫矩阵）result的第一个位置
T=70;
H='bear.jpg';
W='ldu3232.jpg';
%W='200100732.jpg';
Host=imread(H);%读入宿主图片
Water=imread(W);%读入水印
lenw=size(Water,1)*size(Water,2)*8;%水印矩阵行数*列数，32*32*8水印图像所占的bit位吗
% figure(1),subplot(121),imshow(Host),title('Host image');
% subplot(122),imshow(Water),title('watermark image');%显示水印图片和宿主图片
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
levels = 1 ;        % Decomposition level
pfilter = 'maxflat' ;              % Pyramidal filter
dfilter = 'dmaxflat7' ;  
y = nsctdec(Hostlevel, levels,dfilter,pfilter );%非下contourlet变换，提取其低频分量
lowpass=y{1,1};%每一层得低频分量
Waterold=levelwatermark(level,:);%将对应的水印level二进制序列赋给waterold，水印分了三层，一层一层的嵌入
[Hr,Hc]=size(y{1,1});%返回Hostlevel二维矩阵的行列数。
waterposition=1;%嵌入水印位置，初始值为1
for i=1:Hr/blocksize%将载体图像分块，按行可以分成几个。
        if mod(i,2)==1 && mod(level,2)==1%根据层数和行数构造一个网格选块布局
            star=1;
        else
            star=2;
        end
    for j=star:2:Hc/blocksize%先对列进行划分，列按照4*4分块原则可以分多少次，间隔选列，交叉进行，步长为2
        watermark=Waterold(1,waterposition);%将此时该嵌入的水印的个数赋给水印位
        block=lowpass((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);%分块位置
        h=sum((sum(block))');
        if (watermark=='1') %水印嵌入值等于1的时候，量化过程如下。
           if mod(h,T)<0.25*T
               h1=h-mod(h,T)-0.25*T;
           else
                h1=h-mod(h,T)+0.75*T;
           end
        else %水印嵌入值等于0的时侯
           if mod(h,T)>=0.75*T
               h1=h-mod(h,T)+1.25*T;
           else
                h1=h-mod(h,T)+0.25*T;
           end
        end
       while  max(max(block))>255||min(min(block))<0    %设置一个外部循环，一直检测块中元素是否越界，直到符合要求为止
       if max(max(block))>255%判断块的最大值是否越界
           block=lowpass((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);%分块位置 
           if h1>h  
               h1=h1-T;
           end
       end
       if min(min(block))<0
           block=lowpass((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);%分块位置
               h1=h1+T;
       end
      d=(h1-h)/20;
       %再次判断数组是否越界
       for m=1:4
            flag=0;%设置一个标志位，当标志位改变时跳出双重循环
           for n=1:4
               if (m==2&&n==2)||(m==2&&n==3)||(m==3&&n==2)||(m==3&&n==3)
                   block(m,n)=block(m,n)+2*d;
               else
                   block(m,n)=block(m,n)+d;
               end
               if block(m,n)>255||block(m,n)<0
                   flag=1;
                   break;
               end
           end
           if flag==1
               break;
           end   
       end
      end
end
        lowpass((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize)=block;%将分块值重新赋予原来的分块
         waterposition=waterposition+1;%水印嵌入位置（嵌入量）加1
        if waterposition>lenw %%%==2799%%判断是否将所有的水印序列都嵌入进去
            break;
        end
     if waterposition>lenw %%==5778==5578==2799
           break;
     end
end            %这一层水印图像嵌入完成
 y{1}=lowpass;
imrec = nsctrec( y,dfilter, pfilter ) ;
watermarkedim(:,:,level)=imrec;%将重构后的载体图像二维矩阵赋给watermarkedim(:,:,level)，测试攻击的时候就用这个名称
imwrite(uint8(imrec),'www.bmp');%将Hostlevel（嵌入水印后的宿主图像）（原来为double类型的）转换成uint8,保存图片.png格式
imwrite(Host(:,:,level),'hhh.bmp');%保存宿主这一层的图像
psnrval=psnrval+PSNR('www.bmp','hhh.bmp');%将这三层加一块，计算两者图片之间的相似度
 end          %RGB三层图像均已嵌入完成
hfx=watermarkedim;
imwrite(uint8(watermarkedim),'watermarked.bmp');%保存嵌入水印的图片为watermarked.bmp（后期攻击所用图片）
%psnrval=PSNR(H,'watermarked.bmp');
hand=figure(2),subplot(131),imshow(Host),title('Original image');
subplot(132),imshow('watermarked.bmp'),title([' unattack PSNR=',num2str(psnrval/3)]);%显示嵌入完成后的水印图片
toc;
%--------------------------------------------------
%水印提取
watermarkedim=double(imread('watermarked.bmp'));
watermarkedim=hfx;%带水印的图像
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
        h2=sum((sum(block))');
        if mod(h2,T)>0.5*T
            watermark='1';
        else
            watermark='0';
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
% result(ss,2)=psnrval/3;
% result(ss,3)=ncval;
