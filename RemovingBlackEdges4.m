function correctedwatermarkedimg=RemovingBlackEdges4(havingblackwatermarkedimg,M,N)
% 容错率高
%%-------裁剪边缘，去掉黑边并缩放到指定大小[M,N]-------%%
[f,c]=find(havingblackwatermarkedimg(:,:,1));
minf=min(f);maxf=max(f);
minc=min(c);maxc=max(c);

%阈值设置不合理，其值设置应该与相邻行或列的元素有比例关系，例如：首行占下行同列元素的70%~80%方可确定为上边界
deta=30; ratio=0.618; %依据：黄金分割点
hr=havingblackwatermarkedimg(:,:,1);

%测试图片左边非零元素的个数
mincc=minc;
for i=minc:maxc
    if length(find(hr(minf:maxf,i)>=deta))>(maxf-minf+1)*ratio
        mincc=i;
        break;
    end
end

%测试图片右边非零元素的个数
maxcc=maxc;
for i=maxc:-1:minc
    if length(find(hr(minf:maxf,i)>=deta))>(maxf-minf+1)*ratio
        maxcc=i;
        break;
    end
end

%测试图片上边非零元素的个数
minff=minf;
for i=minf:maxf
    if length(find(hr(i,mincc:maxcc)>=deta))>(maxcc-mincc+1)*ratio
        minff=i;
        break;
    end
end

%测试图片下边非零元素的个数
maxff=maxf;
for i=maxf:-1:minf
    if length(find(hr(i,mincc:maxcc)>=deta))>(maxcc-mincc+1)*ratio
        maxff=i;
        break;
    end
end

I2=havingblackwatermarkedimg(minff:maxff,mincc:maxcc,:);
% figure(2),imshow(I2),title('去黑边');
% imshow(uint8(I2)),title('I2');

correctedwatermarkedimg=imresize(I2,[M N]);%%%对裁边的图像进行缩放到规则尺寸
end
