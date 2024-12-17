function correctedwatermarkedimg=RemovingBlackEdges4(havingblackwatermarkedimg,M,N)
% �ݴ��ʸ�
%%-------�ü���Ե��ȥ���ڱ߲����ŵ�ָ����С[M,N]-------%%
[f,c]=find(havingblackwatermarkedimg(:,:,1));
minf=min(f);maxf=max(f);
minc=min(c);maxc=max(c);

%��ֵ���ò�������ֵ����Ӧ���������л��е�Ԫ���б�����ϵ�����磺����ռ����ͬ��Ԫ�ص�70%~80%����ȷ��Ϊ�ϱ߽�
deta=30; ratio=0.618; %���ݣ��ƽ�ָ��
hr=havingblackwatermarkedimg(:,:,1);

%����ͼƬ��߷���Ԫ�صĸ���
mincc=minc;
for i=minc:maxc
    if length(find(hr(minf:maxf,i)>=deta))>(maxf-minf+1)*ratio
        mincc=i;
        break;
    end
end

%����ͼƬ�ұ߷���Ԫ�صĸ���
maxcc=maxc;
for i=maxc:-1:minc
    if length(find(hr(minf:maxf,i)>=deta))>(maxf-minf+1)*ratio
        maxcc=i;
        break;
    end
end

%����ͼƬ�ϱ߷���Ԫ�صĸ���
minff=minf;
for i=minf:maxf
    if length(find(hr(i,mincc:maxcc)>=deta))>(maxcc-mincc+1)*ratio
        minff=i;
        break;
    end
end

%����ͼƬ�±߷���Ԫ�صĸ���
maxff=maxf;
for i=maxf:-1:minf
    if length(find(hr(i,mincc:maxcc)>=deta))>(maxcc-mincc+1)*ratio
        maxff=i;
        break;
    end
end

I2=havingblackwatermarkedimg(minff:maxff,mincc:maxcc,:);
% figure(2),imshow(I2),title('ȥ�ڱ�');
% imshow(uint8(I2)),title('I2');

correctedwatermarkedimg=imresize(I2,[M N]);%%%�Բñߵ�ͼ��������ŵ�����ߴ�
end
