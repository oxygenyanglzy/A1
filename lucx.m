clc;
clear all;
blocksize=4;
% %�����������ң����϶��µ�˳�������ѡ��Ƕ��� ����������DC����������
% T=100;
H='Barbara.jpg';
%W='ldu3232.jpg';
W='200100732.jpg';
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
%Hostlevelback=Hostlevel;
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
        [l,u]=lu(block);
         h=sum(u(1,:));
         if level==1
            T=70;
        else if level==2
                T=70*0.94;
            else 
                T=70*0.78;
            end
        end
       
%        T=70;
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
      d=(h1-h)/4;
           for n=1:4
                   u(1,n)=u(1,n)+d;
           end
       blocknew=l*u;
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
        [l2,u2]=lu(blocknew);
        h2=sum(u2(1,:));
        if level==1
            T=70;
        else if level==2
                T=70*0.94;
            else 
                T=70*0.78;
            end
        end
%        T=70;
        if mod(h2,T)>0.5*T
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