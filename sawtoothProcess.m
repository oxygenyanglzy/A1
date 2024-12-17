function I1=sawtoothProcess(I1)
%先判断个别像素是哪个(对于一行，从左往右相邻两个数挨个比较，如果差值小于diffval，个别像素为右边的像素（一定是右边吗）)，再更改像素
%从左往右相邻两个数挨个比较原因（不用两数相除看相似性原因）：都是差两个像素值134/136=0.9853，4/6=0.6667，结果差距大不能这么比
%只要三个通道 有一个通道相邻两数差距大，就能确定是哪些个别像素位置，然后把三个通道的都改了

%只看一个通道大于阈值不可，有些锯齿只有一个通道的像素不连续，因此只改一个通道的就行了
%有些地方就是颜色变化大的中间界限，不用改，此时跳过，判断条件可以是当前像素和后面五个(若两部分锯齿之间正常像素只有一两个该怎么办)像素相比都大于阈值
%三个通道的阈值是否应该不同？

%有些地方就是颜色变化大的中间界限，所以不能左右两个数相比较了，选择上下两个数进行比较，那么列也上下比较
%第一个和最后一个像素（图像四个顶点像素）未考虑比较后差值较大的情况
[mI1,nI1,~]=size(I1);  sawtoothThres=10;

    
% 上行
i=2;
while i<nI1
    k=1;
    if abs(I1(2,i,1)-I1(1,i,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
        sawtoothPos(k)=i; k=k+1; j=i+1;
        while j<nI1 && abs(I1(2,j,1)-I1(1,j,1))>sawtoothThres  %||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
            sawtoothPos(k)=j; k=k+1; j=j+1;
        end
        k=k-1;%有锯齿效应的像素个数
        %线性插值--消除锯齿效应
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

% 下行
i=2;
while i<nI1
    k=1;
    if abs(I1(mI1-1,i,1)-I1(mI1,i,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
        sawtoothPos(k)=i; k=k+1; j=i+1;
        while j<nI1 && abs(I1(mI1-1,j,1)-I1(mI1,j,1))>sawtoothThres  %||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
            sawtoothPos(k)=j; k=k+1; j=j+1;
        end
        k=k-1;%有锯齿效应的像素个数
        %线性插值--消除锯齿效应
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

%  左列
i=1;
while i<mI1-1 
    if abs(I1(i,1,1)-I1(i+1,1,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
        sawtoothPos(1)=i+1; j=i+1+1;
        while j<mI1-1 && abs(I1(j,1,1)-I1(j+1,1,1))>sawtoothThres  %||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
            j=j+1;
        end
        sawtoothPos(2)=j-1; 
        k=sawtoothPos(2)-sawtoothPos(1)+1;%有锯齿效应的像素个数
         
        %线性插值--消除锯齿效应
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

%  右列
i=1;
while i<mI1-1 
    if abs(I1(i,nI1,1)-I1(i+1,nI1,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
        sawtoothPos(1)=i+1; j=i+1+1;
        while j<mI1-1 && abs(I1(j,nI1,1)-I1(j+1,nI1,1))>sawtoothThres  %||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
            j=j+1;
        end
        sawtoothPos(2)=j-1; 
        k=sawtoothPos(2)-sawtoothPos(1)+1;%有锯齿效应的像素个数
         
        %线性插值--消除锯齿效应
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

%左右比较
% % 左列
% i=1;
% while i<=mI1
%     k=1;
%     if abs(I1(i,1,1)-I1(i,2,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
%         sawtoothPos(k)=i; k=k+1; j=i+1;
%         while abs(I1(j,1,1)-I1(j,2,1))>sawtoothThres && j<=mI1%||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
%             sawtoothPos(k)=j; k=k+1; j=j+1;
%         end
%         k=k-1;%有锯齿效应的像素个数
%         %线性插值--消除锯齿效应
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
% % 右列
% i=2;
% while i<mI1
%     k=1;
%     if abs(I1(i,nI1-1,1)-I1(i,nI1,1))>sawtoothThres%||abs(I1(1,i,2)-I1(1,i+1,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,i+1,3))>sawtoothThres
%         sawtoothPos(k)=i; k=k+1; j=i+1;
%         while abs(I1(j,nI1-1,1)-I1(j,nI1,1))>sawtoothThres && j<=mI1%||abs(I1(1,i,2)-I1(1,j,2))>sawtoothThres||abs(I1(1,i,3)-I1(1,j,3))>sawtoothThres
%             sawtoothPos(k)=j; k=k+1; j=j+1;
%         end
%         k=k-1;%有锯齿效应的像素个数
%         %线性插值--消除锯齿效应
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