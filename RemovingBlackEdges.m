function correctedwatermarkedimg=RemovingBlackEdges(havingblackwatermarkedimg,M,N)
%2019.12.14���  �ϳ�����뷨
%%-------�ü���Ե��ȥ���ڱ߲����ŵ�ָ����С[M,N]-------%%
[f,c]=find(havingblackwatermarkedimg(:,:,1));
minf=min(f);maxf=max(f);
minc=min(c);maxc=max(c);

%����ͼƬ��߷���Ԫ�صĸ���
count=0;

%��ֵ���ò�������ֵ����Ӧ���������л��е�Ԫ���б�����ϵ�����磺����ռ����ͬ��Ԫ�ص�70%~80%����ȷ��Ϊ�ϱ߽�
deta=30;


mincc=minc;
for j=minc:maxc
    for i=minf:maxf
        if havingblackwatermarkedimg(i,j,1)>=deta % && havingblackwatermarkedimg(i,j,2)>=deta && havingblackwatermarkedimg(i,j,1)>=deta
            count=count+1;
        end
    end
    if count<=(maxf-minf)/2  %�������Ԫ�ذ���������<30��Ԫ�أ���˵��������Ч������+1��ֱ����ЧΪֹ
        mincc=mincc+1;
    else
        break;
    end
end
%����ͼƬ�ұ߷���Ԫ�صĸ���
count=0;
maxcc=maxc;
for j=minc:maxc
    jj=maxc+minc-j;
    for i=minf:maxf
        if havingblackwatermarkedimg(i,jj,1)>=deta %&& havingblackwatermarkedimg(i,jj,2)>=deta && havingblackwatermarkedimg(i,jj,1)>=deta
            count=count+1;
        end
    end
    if count<=(maxf-minf)/2 %���β��Ԫ�ذ���������<30��Ԫ�أ���˵��������Ч������-1��ֱ����ЧΪֹ
        maxcc=maxcc-1;
    else
        break;
    end
end
%����ͼƬ�ϱ߷���Ԫ�صĸ���
count=0;
minff=minf;
for i=minf:maxf
    for j=minc:maxc
        if havingblackwatermarkedimg(i,j,1)>=deta% && havingblackwatermarkedimg(i,j,2)>=deta && havingblackwatermarkedimg(i,j,1)>=deta
            count=count+1;
        end
    end
    if count<=(maxc-minc)/2 %�������Ԫ�ذ���������<30��Ԫ�أ���˵��������Ч������+1��ֱ����ЧΪֹ
        minff=minff+1;
    else
        break;
    end
end
%����ͼƬ�±߷���Ԫ�صĸ���
count=0;
maxff=maxf;
for i=minf:maxf
    ii=minf+maxf-i;
    for j=minc:maxc
        if havingblackwatermarkedimg(ii,j,1)>=deta% && havingblackwatermarkedimg(ii,j,2)>=deta && havingblackwatermarkedimg(ii,j,1)>=deta
            count=count+1;
        end
    end
    if count<=(maxc-minc)/2 %���β��Ԫ�ذ���������<30��Ԫ�أ���˵��������Ч������-1��ֱ����ЧΪֹ
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
correctedwatermarkedimg=imresize(I1,[M N]);%%%�Բñߵ�ͼ��������ŵ�����ߴ�
end


%2019.12.09���  ������뷨
%%-------�ü���Ե��ȥ���ڱ߲����ŵ�ָ����С[M,N]-------%%
% [r3,c3,~]=size(havingblackwatermarkedimg);
% flag=0; %��λͼ��������
% for j=1:c3
%     for i=1:r3
%         if havingblackwatermarkedimg(i,j,3)>=30  %%���еķ������ҵ���һ����Ϊ0 �ĵ������[rr1,cc1]
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
% flag=0;%��λͼ��������
% for i=1:r3
%     for j=1:c3
%         if havingblackwatermarkedimg(i,j,3)>=30  %%���еķ������ҵ���һ����Ϊ0 �ĵ������[rr2,cc2]
%             rr2=i; %���з����϶�λ����� �����Ե��
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
% %��Ϊ��ת���ͼ���ǶԳƵģ����Ի������˷�ͼ�񲿷ֵĿ����һ����
% rend=r3-r+1;
% cstart=c;
% cend=c3-c+1;
% %����ͼƬ��߷���Ԫ�صĸ���
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
% %����ͼƬ�ұ߷���Ԫ�صĸ���
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
% %����ͼƬ�ϱ߷���Ԫ�صĸ���
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
% %����ͼƬ�±߷���Ԫ�صĸ���
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
% correctedwatermarkedimg=imresize(I1,[M N]);%%%�Բñߵ�ͼ��������ŵ�����ߴ�

