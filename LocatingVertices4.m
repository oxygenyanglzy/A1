function [cornerpoint,attacktype]=LocatingVertices4(attackedimg)
%像素检测，有局限性，若边缘处有一大块黑的无法检测，边缘就断了
%还有很多情况没考虑，还需优化
ai=attackedimg(:,:,1)+attackedimg(:,:,2)+attackedimg(:,:,3);
ei=logical(ai);
% figure(2),imshow(ei);
% imwrite(ei,'logicalimg.bmp');
% tic
[m,n]=size(ei);
%boder=5; if sum(sum(ei(1:boder,:))) %如果只平移(3,3)，那也检测不出来，只能检测到 测不出边才可以
flag1=0; flag2=0; attacktype=0;
%左边处理检测不到的边
f=find(ai(:,1));
if length(f) > m/2 %如果是小于m/2的图得到的水印也是看不出来了，所以不考虑
    a=f(1); b=f(length(f)); %都这样了，直接给顶点赋值得了呗
    %     ei(a:b,1)=1;
    cornerpoint(:,:,1)=[1 a]; cornerpoint(:,:,4)=[1 b];
    flag1=1;
end
%右边处理检测不到的边
f=find(ai(:,n));
if length(f) > m/2
    a=f(1); b=f(length(f)); %     ei(a:b,1)=1;
    cornerpoint(:,:,2)=[n a]; cornerpoint(:,:,3)=[n b];
    flag2=1;
end
% flag1和flag2都检测不到边，两种可能：水平仿射和平移，还需要再判断
% flag1和flag2只有一个检测不到边，就是平移
% 平移不需要这么找点，横坐标直接赋值，纵坐标只需要判断哪一列有值即可
if flag1 && flag2
    %需要判断是垂直仿射还是平移
    if abs(cornerpoint(1,2,1) - cornerpoint(1,2,2))<=1 %是平移, 平移不适合用透视变换，水平仿射用透视变换
        attacktype=1;
    end
    return;
elseif flag1 || flag2
    %是平移
    attacktype=1;
    if flag1
        c=find(ai(ceil(m/2),:),1,'last');
        cornerpoint(:,:,2)=[c a]; cornerpoint(:,:,3)=[c b];
    else
        c=find(ai(ceil(m/2),:),1);
        cornerpoint(:,:,1)=[c a]; cornerpoint(:,:,4)=[c b];
    end
    return;
end

rstep=2;
%---------左上角---------
rmid=ceil(m/2); %中间行
rmidcol=ei(rmid,:);
h=10;
%左边的两个角（旋转<45°，但是倾斜的情况一般不会旋转这么大角度）
rmidc=find(rmidcol,1); %找到中间行第一个'1'像素
%旋转45°情况
if rmidc==1
    cornerpoint(:,:,1)=[rmidc rmid]; cornerpoint(:,:,2)=[rmid rmidc];
    cornerpoint(:,:,3)=[n rmid]; cornerpoint(:,:,4)=[rmid m];
    return;
end

rmid=rmid-h;

if sum(ei(rmid,1:rmidc-1))==0 %则为右斜
    rmidc1=rmidc-1+find(ei(rmid,rmidc:n),1);
    diff=abs(rmidc1-rmidc);
    
    while rmidc1<=n && rmid>=1
        rmid=rmid-h; rmidcol=rmidc1+diff-rstep:rmidc1+diff+rstep;
        if rmid<1
            rmid=rmid+h; diff=ceil(diff/h); h=1;
        elseif sum(ei(rmid,rmidcol))==0
            if rmid-1>0 && sum(ei(rmid-1,rmidcol))
                rmid=rmid-1;
                rmidc1=rmidcol(find(ei(rmid,rmidcol),1));
            elseif h==1
                rmid=rmid+h;
                col=rmidcol(find(ei(rmid,rmidcol),1));
                cornerpoint(:,:,1)=[col,rmid];
                break;
            else
                rmid=rmid+h; diff=ceil(diff/h); h=1;
            end
        elseif rmid==1
            col=rmidcol(find(ei(rmid,rmidcol),1));
            cornerpoint(:,:,1)=[col,rmid];
            break;
        else
            rmidc1=rmidcol(find(ei(rmid,rmidcol),1));
        end
    end
    
else %左斜
    if sum(ei(rmid-h,1:rmidc))==0 %旋转44°会遇到的情况
        rmid=rmid+h-1; h=1;
    end
    
    rmidc1=find(ei(rmid,1:rmidc),1); % 新行rmid，新列rmidc1(这个肯定不会超)
    diff=abs(rmidc1-rmidc);
    %     h=2*h; diff=2*diff; % h=10 （tanθ=diff/h=2diff/2h=(diff/h)/1）
    while rmidc1>=1 && rmid>=1 %还没考虑超出四边的框的情况，左边框考虑了，上边框只要超过了回溯即可，右和下同理
        % 1如果垂直移动h像素后判断水平范围内无像素，则回溯一次，令h=1，如果这移动后水平范围内还是无像素，则再次回溯，直到水平范围内有像素
        % 2最后几步只垂直移动一个像素(h=1)，水平移动diff/h个范围内的像素(diff/h=4/5=0.8∈[0,1]，即0和1个像素)，并判断右边(包括右上、右下)是否有像素，如果有则回溯一次，便是顶点
        % 3如果最后水平范围内有超出边框的点，但是没超出的点的上边和左上还有像素，则令h=1；否则那个没超出的点就是顶点；
        % 4如果最后水平范围内有超出边框的点，但是其他也无像素，则回溯一次；
        % 5如果最后水平范围内所有点超出边框，回溯一次
        rmid=rmid-h; rmidcol=rmidc1-diff-rstep:rmidc1-diff+rstep; %可能第一次上移后的左移就出左边框了
        if rmid<1
            rmid=rmid+h-1; diff=ceil(diff/h); h=1;
        end
        
        if rmidcol(length(rmidcol))<=0 % 5
            %回溯一次，令h=1
            rmid=rmid+h; diff=ceil(diff/h); h=1;
        elseif rmidcol(1)<=0
            if sum(ei(rmid,1:rmidcol(length(rmidcol))))==0 % 4
                if rmid-1>0 && sum(ei(rmid-1,1:rmidcol(length(rmidcol))))
                    rmid=rmid-1;
                    rmidc1=rmidcol(find(ei(rmid,1:rmidcol(length(rmidcol))),1));
                else
                    %回溯一次，令h=1
                    rmid=rmid+h; diff=ceil(diff/h); h=1;
                end
            else %3
                col=find(ei(rmid,1:rmidcol(length(rmidcol))),1);
                while ((rmid-1)>0 && ei(rmid-1,col)) || ((col-1)>0 && ei(rmid-1,col-1))
                    if ((col-1)>0 && ei(rmid-1,col-1))
                        rmid=rmid-1; col=col-1;
                    elseif ((rmid-1)>0 && ei(rmid-1,col))
                        rmid=rmid-1;
                    end
                end
                cornerpoint(:,:,1)=[col,rmid];
                break;
            end
        elseif sum(ei(rmid,rmidcol))==0 % 1
            if (rmid+1)>0 && sum(ei(rmid+1,rmidcol))
                rmid=rmid+1;
            else
                %回溯一次，令h=1
                rmid=rmid+h; diff=ceil(diff/h); h=1;
            end
        elseif rmid==1
            col=rmidcol(find(ei(rmid,rmidcol),1));
            cornerpoint(:,:,1)=[col,rmid];
            break;
        else %正常
            rmidc1=rmidcol(find(ei(rmid,rmidcol),1));
        end
    end
end

%---------左下角---------
rmid=ceil(m/2); %中间行
rmidcol=ei(rmid,:);
h=10;
rmidc=find(rmidcol,1); %找到中间行第一个'1'像素
rmid=rmid+h;
if sum(ei(rmid,1:rmidc-1))==0 %则为右斜
    rmidc1=rmidc-1+find(ei(rmid,rmidc:n),1);
    diff=abs(rmidc1-rmidc);
    
    while rmidc1<=n && rmid<=m
        rmid=rmid+h; rmidcol=rmidc1+diff-rstep:rmidc1+diff+rstep;
        
        if rmid>m
            rmid=rmid-h; diff=ceil(diff/h); h=1;
        elseif sum(ei(rmid,rmidcol))==0
            if rmid+1<=m && sum(ei(rmid+1,rmidcol))
                rmid=rmid+1;
                rmidc1=rmidcol(find(ei(rmid,rmidcol),1));
            elseif h==1
                rmid=rmid-h;
                col=rmidcol(find(ei(rmid,rmidcol),1));
                cornerpoint(:,:,4)=[col,rmid];
                break;
            else
                rmid=rmid-h; diff=ceil(diff/h); h=1;
            end
        elseif rmid==m
            col=rmidcol(find(ei(rmid,rmidcol),1));
            cornerpoint(:,:,4)=[col,rmid];
            break;
        else
            rmidc1=rmidcol(find(ei(rmid,rmidcol),1));
        end
    end
    
else %左斜
    if sum(ei(rmid+h,1:rmidc))==0 %旋转-44°会遇到的情况
        rmid=rmid-h+1; h=1;
    end
    rmidc1=find(ei(rmid,1:rmidc),1);
    diff=abs(rmidc1-rmidc);
    while rmidc1>=1 && rmid<=m
        rmid=rmid+h; rmidcol=rmidc1-diff-rstep:rmidc1-diff+rstep;
        if rmid>m
            rmid=rmid-h+1; diff=ceil(diff/h); h=1;
        end
        if rmidcol(length(rmidcol))<=0 % 5
            %回溯一次，令h=1
            rmid=rmid-h; diff=ceil(diff/h); h=1;
        elseif rmidcol(1)<=0
            if sum(ei(rmid,1:rmidcol(length(rmidcol))))==0 % 4
                if rmid+1<=m && sum(ei(rmid+1,1:rmidcol(length(rmidcol))))
                    rmid=rmid+1;
                    rmidc1=rmidcol(find(ei(rmid,1:rmidcol(length(rmidcol))),1));
                else
                    %回溯一次，令h=1
                    rmid=rmid-h; diff=ceil(diff/h); h=1;
                end
            else %3
                col=find(ei(rmid,1:rmidcol(length(rmidcol))),1);
                while ((rmid+1)<=m && ei(rmid+1,col)) || ((col-1)>0 && ei(rmid+1,col-1))
                    if ((col-1)>0 && ei(rmid+1,col-1))
                        rmid=rmid+1; col=col-1;
                    elseif ((rmid+1)<=m && ei(rmid+1,col))
                        rmid=rmid+1;
                    end
                end
                cornerpoint(:,:,4)=[col,rmid];
                break;
            end
        elseif sum(ei(rmid,rmidcol))==0 % 1
            if (rmid+1)>0 && sum(ei(rmid+1,rmidcol))
                rmid=rmid+1;
            else
                %回溯一次，令h=1
                rmid=rmid-h; diff=ceil(diff/h); h=1;
            end
        elseif rmid==m
            col=rmidcol(find(ei(rmid,rmidcol),1));
            cornerpoint(:,:,4)=[col,rmid];
            break;
        else %正常
            rmidc1=rmidcol(find(ei(rmid,rmidcol),1));
        end
    end
end


% %---------右上角---------
rmid=ceil(m/2); %中间行
rmidcol=ei(rmid,:);
h=10;
rmidc=find(rmidcol,1,'last');
rmid=rmid-h;
if sum(ei(rmid,rmidc+1:n))==0 %则为左斜
    rmidc1=find(ei(rmid,1:rmidc),1,'last');
    diff=abs(rmidc1-rmidc);
    
    while rmidc1<=n && rmid>=1
        rmid=rmid-h; rmidcol=rmidc1-diff-rstep:rmidc1-diff+rstep;
        if rmid<1
            rmid=rmid+h; diff=ceil(diff/h); h=1;
        elseif sum(ei(rmid,rmidcol))==0
            if rmid-1>0 && sum(ei(rmid-1,rmidcol))
                rmid=rmid-1;
                rmidc1=rmidcol(find(ei(rmid,rmidcol),1,'last'));
            elseif h==1
                rmid=rmid+h;
                col=rmidcol(find(ei(rmid,rmidcol),1,'last'));
                cornerpoint(:,:,2)=[col,rmid];
                break;
            else
                rmid=rmid+h; diff=ceil(diff/h); h=1;
            end
        elseif rmid==1
            col=rmidcol(find(ei(rmid,rmidcol),1));
            cornerpoint(:,:,2)=[col,rmid];
            break;
        else
            rmidc1=rmidcol(find(ei(rmid,rmidcol),1,'last'));
        end
    end
    
else %右斜
    if sum(ei(rmid-h,rmidc:n))==0 %旋转-44°会遇到的情况
        rmid=rmid+h-1; h=1;
    end
    rmidc1=rmidc-1+find(ei(rmid,rmidc:n),1,'last');
    diff=abs(rmidc1-rmidc);
    while rmidc1<=n && rmid>=1
        rmid=rmid-h; rmidcol=rmidc1+diff-rstep:rmidc1+diff+rstep;
        if rmid<1
            rmid=rmid+h-1; diff=ceil(diff/h); h=1;
        end
        
        if rmidcol(1)>n
            %回溯一次，令h=1
            rmid=rmid+h; diff=ceil(diff/h); h=1;
        elseif rmidcol(length(rmidcol))>n
            if sum(ei(rmid,rmidcol(1):n))==0
                if rmid-1>0 && sum(ei(rmid-1,rmidcol(1):n))
                    rmid=rmid-1;
                    rmidc1=rmidcol(find(ei(rmid,rmidcol(1):n),1,'last'));
                else
                    %回溯一次，令h=1
                    rmid=rmid+h; diff=ceil(diff/h); h=1;
                end
            else
                col=rmidcol(find(ei(rmid,rmidcol(1):n),1,'last'));
                while ((rmid-1)>0 && ei(rmid-1,col)) || ((col+1)<=n && ei(rmid-1,col+1))
                    if ((col+1)<=n && ei(rmid-1,col+1))
                        rmid=rmid-1; col=col+1;
                    elseif ((rmid-1)>0 && ei(rmid-1,col))
                        rmid=rmid-1;
                    end
                end
                cornerpoint(:,:,2)=[col,rmid];
                break;
            end
        elseif sum(ei(rmid,rmidcol))==0
            if (rmid-1)>0 && sum(ei(rmid-1,rmidcol))
                rmid=rmid-1;
            else
                %回溯一次，令h=1
                rmid=rmid+h; diff=ceil(diff/h); h=1;
            end
        elseif rmid==1
            col=rmidcol(find(ei(rmid,rmidcol),1));
            cornerpoint(:,:,2)=[col,rmid];
            break;
        else
            rmidc1=rmidcol(find(ei(rmid,rmidcol),1,'last'));
        end
    end
end

% %---------右下角---------
rmid=ceil(m/2); %中间行
rmidcol=ei(rmid,:);
h=10;
rmidc=find(rmidcol,1,'last');
rmid=rmid+h;
if sum(ei(rmid,rmidc+1:n))==0 %则为左斜
    rmidc1=find(ei(rmid,1:rmidc),1,'last');
    diff=abs(rmidc1-rmidc);
    
    while rmidc1<=n && rmid>=1
        rmid=rmid+h; rmidcol=rmidc1-diff-rstep:rmidc1-diff+rstep;
        if rmid>m
            rmid=rmid-h; diff=ceil(diff/h); h=1;
        elseif sum(ei(rmid,rmidcol))==0
            if rmid+1<=m && sum(ei(rmid+1,rmidcol))
                rmid=rmid+1;
                rmidc1=rmidcol(find(ei(rmid,rmidcol),1,'last'));
            elseif h==1
                rmid=rmid-h;
                col=rmidcol(find(ei(rmid,rmidcol),1,'last'));
                cornerpoint(:,:,3)=[col,rmid];
                break;
            else
                rmid=rmid-h; diff=ceil(diff/h); h=1;
            end
        elseif rmid==m
            col=rmidcol(find(ei(rmid,rmidcol),1,'last'));
            cornerpoint(:,:,3)=[col,rmid];
            break;
        else
            rmidc1=rmidcol(find(ei(rmid,rmidcol),1,'last'));
        end
    end
    
else %右斜
    if sum(ei(rmid+h,rmidc:n))==0 %旋转44°会遇到的情况
        rmid=rmid-h+1; h=1;
    end
    rmidc1=rmidc-1+find(ei(rmid,rmidc:n),1,'last');
    diff=abs(rmidc1-rmidc);
    while rmidc1<=n && rmid>=1
        rmid=rmid+h; rmidcol=rmidc1+diff-rstep:rmidc1+diff+rstep;
        if rmid>m
            rmid=rmid-h+1; diff=ceil(diff/h); h=1;
        end
        if rmidcol(1)>n
            %回溯一次，令h=1
            rmid=rmid-h; diff=ceil(diff/h); h=1;
        elseif rmidcol(length(rmidcol))>n
            if sum(ei(rmid,rmidcol(1):n))==0
                if rmid+1<=m && sum(ei(rmid+1,rmidcol(1):n))
                    rmid=rmid+1;
                    rmidc1=rmidcol(find(ei(rmid,rmidcol(1):n),1,'last'));
                else
                    %回溯一次，令h=1
                    rmid=rmid-h; diff=ceil(diff/h); h=1;
                end
            else
                col=rmidcol(find(ei(rmid,rmidcol(1):n),1,'last'));
                while ((rmid+1)<=m && ei(rmid+1,col)) || ((col+1)<=n && ei(rmid+1,col+1))
                    if ((col+1)<=n && ei(rmid+1,col+1))
                        rmid=rmid+1; col=col+1;
                    elseif ((rmid+1)<=m && ei(rmid+1,col))
                        rmid=rmid+1;
                    end
                end
                cornerpoint(:,:,3)=[col,rmid];
                break;
            end
        elseif sum(ei(rmid,rmidcol))==0
            if (rmid+1)<n && sum(ei(rmid+1,rmidcol))
                rmid=rmid+1;
            else
                %回溯一次，令h=1
                rmid=rmid-h; diff=ceil(diff/h); h=1;
            end
        elseif rmid==m
            col=rmidcol(find(ei(rmid,rmidcol),1,'last'));
            cornerpoint(:,:,3)=[col,rmid];
            break;
        else
            rmidc1=rmidcol(find(ei(rmid,rmidcol),1,'last'));
        end
    end
end
% Cornerpoint(:,:,1)=cornerpoint(:,:,2);
% Cornerpoint(:,:,2)=cornerpoint(:,:,3);
% Cornerpoint(:,:,3)=cornerpoint(:,:,4);
% Cornerpoint(:,:,4)=cornerpoint(:,:,1);

% toc
