function I1=sawtoothProcess(I1)
%���жϸ����������ĸ�(����һ�У������������������������Ƚϣ������ֵС��diffval����������Ϊ�ұߵ����أ�һ�����ұ���)���ٸ�������
%�����������������������Ƚ�ԭ�򣨲������������������ԭ�򣩣����ǲ���������ֵ134/136=0.9853��4/6=0.6667�������������ô��
%ֻҪ����ͨ�� ��һ��ͨ�������������󣬾���ȷ������Щ��������λ�ã�Ȼ�������ͨ���Ķ�����

%ֻ��һ��ͨ��������ֵ���ɣ���Щ���ֻ��һ��ͨ�������ز����������ֻ��һ��ͨ���ľ�����
%��Щ�ط�������ɫ�仯����м���ޣ����øģ���ʱ�������ж����������ǵ�ǰ���غͺ������(�������־��֮����������ֻ��һ��������ô��)������ȶ�������ֵ
%����ͨ������ֵ�Ƿ�Ӧ�ò�ͬ��

%��Щ�ط�������ɫ�仯����м���ޣ����Բ���������������Ƚ��ˣ�ѡ���������������бȽϣ���ô��Ҳ���±Ƚ�
%��һ�������һ�����أ�ͼ���ĸ��������أ�δ���ǱȽϺ��ֵ�ϴ�����
[mI1,nI1,~]=size(I1);  sawtoothThres=10;

    
% ����
i=2;
while i<nI1
    k=1;
    if abs(I1(2,i,1)-I1(1,i,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
        sawtoothPos(k)=i; k=k+1; j=i+1;
        while j<nI1 && abs(I1(2,j,1)-I1(1,j,1))>sawtoothThres  %||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
            sawtoothPos(k)=j; k=k+1; j=j+1;
        end
        k=k-1;%�о��ЧӦ�����ظ���
        %���Բ�ֵ--�������ЧӦ
        x=[1 k+2];
        xi=2:k+1;
        for lev=1:3
            y=[I1(1,sawtoothPos(1)-1,lev) I1(1,sawtoothPos(k)+1,lev)];
            y1=interp1(x,y,xi);
            I1(1,sawtoothPos(1):sawtoothPos(k),lev)=y1;
        end
%           break;
        i=j+1;
    else
        i=i+1;
    end
end

% ����
i=2;
while i<nI1
    k=1;
    if abs(I1(mI1-1,i,1)-I1(mI1,i,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
        sawtoothPos(k)=i; k=k+1; j=i+1;
        while j<nI1 && abs(I1(mI1-1,j,1)-I1(mI1,j,1))>sawtoothThres  %||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
            sawtoothPos(k)=j; k=k+1; j=j+1;
        end
        k=k-1;%�о��ЧӦ�����ظ���
        %���Բ�ֵ--�������ЧӦ
        x=[1 k+2];
        xi=2:k+1;
        for lev=1:3
            y=[I1(mI1,sawtoothPos(1)-1,lev) I1(mI1,sawtoothPos(k)+1,lev)];
            y1=interp1(x,y,xi);
            I1(mI1,sawtoothPos(1):sawtoothPos(k),lev)=y1;
        end
        %         break;
        i=j+1;
    else
        i=i+1;
    end
end

%  ����
i=1;
while i<mI1-1 
    if abs(I1(i,1,1)-I1(i+1,1,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
        sawtoothPos(1)=i+1; j=i+1+1;
        while j<mI1-1 && abs(I1(j,1,1)-I1(j+1,1,1))>sawtoothThres  %||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
            j=j+1;
        end
        sawtoothPos(2)=j-1; 
        k=sawtoothPos(2)-sawtoothPos(1)+1;%�о��ЧӦ�����ظ���
         
        %���Բ�ֵ--�������ЧӦ
        x=[1 k+2];
        xi=2:k+1;
        for lev=1:3
            y=[I1(sawtoothPos(1)-1,1,lev) I1(sawtoothPos(2)+1,1,lev)];
            y1=interp1(x,y,xi);
            I1(sawtoothPos(1):sawtoothPos(2),1,lev)=y1;
        end
%                 break;
        i=j+1;
    else
        i=i+1;
    end
end

%  ����
i=1;
while i<mI1-1 
    if abs(I1(i,nI1,1)-I1(i+1,nI1,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
        sawtoothPos(1)=i+1; j=i+1+1;
        while j<mI1-1 && abs(I1(j,nI1,1)-I1(j+1,nI1,1))>sawtoothThres  %||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
            j=j+1;
        end
        sawtoothPos(2)=j-1; 
        k=sawtoothPos(2)-sawtoothPos(1)+1;%�о��ЧӦ�����ظ���
         
        %���Բ�ֵ--�������ЧӦ
        x=[1 k+2];
        xi=2:k+1;
        for lev=1:3
            y=[I1(sawtoothPos(1)-1,nI1,lev) I1(sawtoothPos(2)+1,nI1,lev)];
            y1=interp1(x,y,xi);
            I1(sawtoothPos(1):sawtoothPos(2),nI1,lev)=y1;
        end
%                 break;
        i=j+1;
    else
        i=i+1;
    end
end

%���ұȽ�
% % ����
% i=1;
% while i<=mI1
%     k=1;
%     if abs(I1(i,1,1)-I1(i,2,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
%         sawtoothPos(k)=i; k=k+1; j=i+1;
%         while abs(I1(j,1,1)-I1(j,2,1))>sawtoothThres && j<=mI1%||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
%             sawtoothPos(k)=j; k=k+1; j=j+1;
%         end
%         k=k-1;%�о��ЧӦ�����ظ���
%         %���Բ�ֵ--�������ЧӦ
%         x=[1 k+2];
%         xi=2:k+1;
%         for lev=1:3
%             y=[I1(sawtoothPos(1)-1,1,lev) I1(sawtoothPos(k)+1,1,lev)];
%             y1=interp1(x,y,xi);
%             I1(sawtoothPos(1):sawtoothPos(k),1,lev)=y1;
%         end
%         %         break;
%         i=j+1;
%     else
%         i=i+1;
%     end
% end
% 
% % ����
% i=2;
% while i<mI1
%     k=1;
%     if abs(I1(i,nI1-1,1)-I1(i,nI1,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
%         sawtoothPos(k)=i; k=k+1; j=i+1;
%         while abs(I1(j,nI1-1,1)-I1(j,nI1,1))>sawtoothThres && j<=mI1%||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
%             sawtoothPos(k)=j; k=k+1; j=j+1;
%         end
%         k=k-1;%�о��ЧӦ�����ظ���
%         %���Բ�ֵ--�������ЧӦ
%         x=[1 k+2];
%         xi=2:k+1;
%         for lev=1:3
%             y=[I1(sawtoothPos(1)-1,nI1,lev) I1(sawtoothPos(k)+1,nI1,lev)];
%             y1=interp1(x,y,xi);
%             I1(sawtoothPos(1):sawtoothPos(k),nI1,lev)=y1;
%         end
%         %         break;
%         i=j+1;
%     else
%         i=i+1;
%     end
% end