function correctedwatermarkedimg=RemovingBlackEdges(havingblackwatermarkedimg,M,N)
%2019.12.14完成  较成熟的想法
%%-------裁剪边缘，去掉黑边并缩放到指定大小[M,N]-------%%
[f,c]=find(havingblackwatermarkedimg(:,:,1));
minf=min(f);maxf=max(f);
minc=min(c);maxc=max(c);

%测试图片左边非零元素的个数
count=0;

%阈值设置不合理，其值设置应该与相邻行或列的元素有比例关系，例如：首行占下行同列元素的70%~80%方可确定为上边界
deta=30;


mincc=minc;
for j=minc:maxc
    for i=minf:maxf
        if havingblackwatermarkedimg(i,j,1)>=deta % && havingblackwatermarkedimg(i,j,2)>=deta && havingblackwatermarkedimg(i,j,1)>=deta
            count=count+1;
        end
    end
    if count<=(maxf-minf)/2  %如果首列元素半数以上是<30的元素，则说明该列无效，列数+1，直至有效为止
        mincc=mincc+1;
    else
        break;
    end
end
%测试图片右边非零元素的个数
count=0;
maxcc=maxc;
for j=minc:maxc
    jj=maxc+minc-j;
    for i=minf:maxf
        if havingblackwatermarkedimg(i,jj,1)>=deta %&& havingblackwatermarkedimg(i,jj,2)>=deta && havingblackwatermarkedimg(i,jj,1)>=deta
            count=count+1;
        end
    end
    if count<=(maxf-minf)/2 %如果尾列元素半数以上是<30的元素，则说明该列无效，列数-1，直至有效为止
        maxcc=maxcc-1;
    else
        break;
    end
end
%测试图片上边非零元素的个数
count=0;
minff=minf;
for i=minf:maxf
    for j=minc:maxc
        if havingblackwatermarkedimg(i,j,1)>=deta% && havingblackwatermarkedimg(i,j,2)>=deta && havingblackwatermarkedimg(i,j,1)>=deta
            count=count+1;
        end
    end
    if count<=(maxc-minc)/2 %如果首行元素半数以上是<30的元素，则说明该行无效，行数+1，直至有效为止
        minff=minff+1;
    else
        break;
    end
end
%测试图片下边非零元素的个数
count=0;
maxff=maxf;
for i=minf:maxf
    ii=minf+maxf-i;
    for j=minc:maxc
        if havingblackwatermarkedimg(ii,j,1)>=deta% && havingblackwatermarkedimg(ii,j,2)>=deta && havingblackwatermarkedimg(ii,j,1)>=deta
            count=count+1;
        end
    end
    if count<=(maxc-minc)/2 %如果尾行元素半数以上是<30的元素，则说明该行无效，行数-1，直至有效为止
        maxff=maxff-1;
    else
        break;
    end
end
I1=havingblackwatermarkedimg(minff:maxff,mincc:maxcc,:);
% I8=havingblackwatermarkedimg(59:571,59:569,:);
% I2=I1(:,:,1);
% I3=I1(:,:,2);
% I4=I1(:,:,3);
% I5=havingblackwatermarkedimg(:,:,1);
% I6=havingblackwatermarkedimg(:,:,2);
% I7=havingblackwatermarkedimg(:,:,3);
correctedwatermarkedimg=imresize(I1,[M N]);%%%对裁边的图像进行缩放到规则尺寸
end


%2019.12.09完成  较早的想法
%%-------裁剪边缘，去掉黑边并缩放到指定大小[M,N]-------%%
% [r3,c3,~]=size(havingblackwatermarkedimg);
% flag=0; %定位图像最左列
% for j=1:c3
%     for i=1:r3
%         if havingblackwatermarkedimg(i,j,3)>=30  %%在列的方向上找到第一个不为0 的点的坐标[rr1,cc1]
%             cc1=j;
%             flag=1;
%             break;
%         end
%     end
%     if flag==1
%         break;
%     end
% end
% 
% flag=0;%定位图像最上行
% for i=1:r3
%     for j=1:c3
%         if havingblackwatermarkedimg(i,j,3)>=30  %%在行的方向上找到第一个不为0 的点的坐标[rr2,cc2]
%             rr2=i; %在列方向上定位最左边 （左边缘）
%             flag=1;
%             break;
%         end
%     end
%     if flag==1
%         break;
%     end
% end
% r=rr2;
% c=cc1;
% 
% rstart=r;
% %因为旋转后的图像是对称的，所以画布两端非图像部分的宽度是一样的
% rend=r3-r+1;
% cstart=c;
% cend=c3-c+1;
% %测试图片左边非零元素的个数
% count=0;
% for j=cstart:cend
%     for i=rstart:rend
%         if havingblackwatermarkedimg(i,j,3)>=30 && havingblackwatermarkedimg(i,j,2)>=30 && havingblackwatermarkedimg(i,j,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(rend-rstart)/2
%         cstart=cstart+1;
%     else
%         break;
%     end
% end
% %测试图片右边非零元素的个数
% count=0;
% for j=cstart:cend
%     jj=cend+cstart-j;
%     for i=rstart:rend
%         if havingblackwatermarkedimg(i,jj,3)>=30 && havingblackwatermarkedimg(i,jj,2)>=30 && havingblackwatermarkedimg(i,jj,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(rend-rstart)/2
%         cend=cend-1;
%     else
%         break;
%     end
% end
% %测试图片上边非零元素的个数
% count=0;
% for i=rstart:rend
%     for j=cstart:cend
%         if havingblackwatermarkedimg(i,j,3)>=30 && havingblackwatermarkedimg(i,j,2)>=30 && havingblackwatermarkedimg(i,j,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(cend-cstart)/2
%         rstart=rstart+1;
%     else
%         break;
%     end
% end
% %测试图片下边非零元素的个数
% count=0;
% for i=rstart:rend
%     ii=rstart+rend-i;
%     for j=cstart:cend
%         if havingblackwatermarkedimg(ii,j,3)>=30 && havingblackwatermarkedimg(ii,j,2)>=30 && havingblackwatermarkedimg(ii,j,1)>=30
%             count=count+1;
%         end
%     end
%     if count<=(cend-cstart)/2
%         rend=rend-1;
%     else
%         break;
%     end
% end
% 
% I1=havingblackwatermarkedimg(rstart:rend,cstart:cend,:);
% correctedwatermarkedimg=imresize(I1,[M N]);%%%对裁边的图像进行缩放到规则尺寸

