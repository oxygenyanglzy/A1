function [cornerpoint,attacktype]=LocatingVertices4(attackedimg)
%���ؼ�⣬�о����ԣ�����Ե����һ���ڵ��޷���⣬��Ե�Ͷ���
%���кܶ����û���ǣ������Ż�
ai=attackedimg(:,:,1)+attackedimg(:,:,2)+attackedimg(:,:,3);
ei=logical(ai);
% figure(2),imshow(ei);
% imwrite(ei,'logicalimg.bmp');
% tic
[m,n]=size(ei);
%boder=5; if sum(sum(ei(1:boder,:))) %���ֻƽ��(3,3)����Ҳ��ⲻ������ֻ�ܼ�⵽ �ⲻ���߲ſ���
flag1=0; flag2=0; attacktype=0;
%��ߴ����ⲻ���ı�
f=find(ai(:,1));
if length(f) > m/2 %�����С��m/2��ͼ�õ���ˮӡҲ�ǿ��������ˣ����Բ�����
    a=f(1); b=f(length(f)); %�������ˣ�ֱ�Ӹ����㸳ֵ������
    %     ei(a:b,1)=1;
    cornerpoint(:,:,1)=[1 a]; cornerpoint(:,:,4)=[1 b];
    flag1=1;
end
%�ұߴ����ⲻ���ı�
f=find(ai(:,n));
if length(f) > m/2
    a=f(1); b=f(length(f)); %     ei(a:b,1)=1;
    cornerpoint(:,:,2)=[n a]; cornerpoint(:,:,3)=[n b];
    flag2=1;
end
% flag1��flag2����ⲻ���ߣ����ֿ��ܣ�ˮƽ�����ƽ�ƣ�����Ҫ���ж�
% flag1��flag2ֻ��һ����ⲻ���ߣ�����ƽ��
% ƽ�Ʋ���Ҫ��ô�ҵ㣬������ֱ�Ӹ�ֵ��������ֻ��Ҫ�ж���һ����ֵ����
if flag1 && flag2
    %��Ҫ�ж��Ǵ�ֱ���仹��ƽ��
    if abs(cornerpoint(1,2,1) - cornerpoint(1,2,2))<=1 %��ƽ��, ƽ�Ʋ��ʺ���͸�ӱ任��ˮƽ������͸�ӱ任
        attacktype=1;
    end
    return;
elseif flag1 || flag2
    %��ƽ��
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
%---------���Ͻ�---------
rmid=ceil(m/2); %�м���
rmidcol=ei(rmid,:);
h=10;
%��ߵ������ǣ���ת<45�㣬������б�����һ�㲻����ת��ô��Ƕȣ�
rmidc=find(rmidcol,1); %�ҵ��м��е�һ��'1'����
%��ת45�����
if rmidc==1
    cornerpoint(:,:,1)=[rmidc rmid]; cornerpoint(:,:,2)=[rmid rmidc];
    cornerpoint(:,:,3)=[n rmid]; cornerpoint(:,:,4)=[rmid m];
    return;
end

rmid=rmid-h;

if sum(ei(rmid,1:rmidc-1))==0 %��Ϊ��б
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
    
else %��б
    if sum(ei(rmid-h,1:rmidc))==0 %��ת44������������
        rmid=rmid+h-1; h=1;
    end
    
    rmidc1=find(ei(rmid,1:rmidc),1); % ����rmid������rmidc1(����϶����ᳬ)
    diff=abs(rmidc1-rmidc);
    %     h=2*h; diff=2*diff; % h=10 ��tan��=diff/h=2diff/2h=(diff/h)/1��
    while rmidc1>=1 && rmid>=1 %��û���ǳ����ıߵĿ���������߿����ˣ��ϱ߿�ֻҪ�����˻��ݼ��ɣ��Һ���ͬ��
        % 1�����ֱ�ƶ�h���غ��ж�ˮƽ��Χ�������أ������һ�Σ���h=1��������ƶ���ˮƽ��Χ�ڻ��������أ����ٴλ��ݣ�ֱ��ˮƽ��Χ��������
        % 2��󼸲�ֻ��ֱ�ƶ�һ������(h=1)��ˮƽ�ƶ�diff/h����Χ�ڵ�����(diff/h=4/5=0.8��[0,1]����0��1������)�����ж��ұ�(�������ϡ�����)�Ƿ������أ�����������һ�Σ����Ƕ���
        % 3������ˮƽ��Χ���г����߿�ĵ㣬����û�����ĵ���ϱߺ����ϻ������أ�����h=1�������Ǹ�û�����ĵ���Ƕ��㣻
        % 4������ˮƽ��Χ���г����߿�ĵ㣬��������Ҳ�����أ������һ�Σ�
        % 5������ˮƽ��Χ�����е㳬���߿򣬻���һ��
        rmid=rmid-h; rmidcol=rmidc1-diff-rstep:rmidc1-diff+rstep; %���ܵ�һ�����ƺ�����ƾͳ���߿���
        if rmid<1
            rmid=rmid+h-1; diff=ceil(diff/h); h=1;
        end
        
        if rmidcol(length(rmidcol))<=0 % 5
            %����һ�Σ���h=1
            rmid=rmid+h; diff=ceil(diff/h); h=1;
        elseif rmidcol(1)<=0
            if sum(ei(rmid,1:rmidcol(length(rmidcol))))==0 % 4
                if rmid-1>0 && sum(ei(rmid-1,1:rmidcol(length(rmidcol))))
                    rmid=rmid-1;
                    rmidc1=rmidcol(find(ei(rmid,1:rmidcol(length(rmidcol))),1));
                else
                    %����һ�Σ���h=1
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
                %����һ�Σ���h=1
                rmid=rmid+h; diff=ceil(diff/h); h=1;
            end
        elseif rmid==1
            col=rmidcol(find(ei(rmid,rmidcol),1));
            cornerpoint(:,:,1)=[col,rmid];
            break;
        else %����
            rmidc1=rmidcol(find(ei(rmid,rmidcol),1));
        end
    end
end

%---------���½�---------
rmid=ceil(m/2); %�м���
rmidcol=ei(rmid,:);
h=10;
rmidc=find(rmidcol,1); %�ҵ��м��е�һ��'1'����
rmid=rmid+h;
if sum(ei(rmid,1:rmidc-1))==0 %��Ϊ��б
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
    
else %��б
    if sum(ei(rmid+h,1:rmidc))==0 %��ת-44������������
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
            %����һ�Σ���h=1
            rmid=rmid-h; diff=ceil(diff/h); h=1;
        elseif rmidcol(1)<=0
            if sum(ei(rmid,1:rmidcol(length(rmidcol))))==0 % 4
                if rmid+1<=m && sum(ei(rmid+1,1:rmidcol(length(rmidcol))))
                    rmid=rmid+1;
                    rmidc1=rmidcol(find(ei(rmid,1:rmidcol(length(rmidcol))),1));
                else
                    %����һ�Σ���h=1
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
                %����һ�Σ���h=1
                rmid=rmid-h; diff=ceil(diff/h); h=1;
            end
        elseif rmid==m
            col=rmidcol(find(ei(rmid,rmidcol),1));
            cornerpoint(:,:,4)=[col,rmid];
            break;
        else %����
            rmidc1=rmidcol(find(ei(rmid,rmidcol),1));
        end
    end
end


% %---------���Ͻ�---------
rmid=ceil(m/2); %�м���
rmidcol=ei(rmid,:);
h=10;
rmidc=find(rmidcol,1,'last');
rmid=rmid-h;
if sum(ei(rmid,rmidc+1:n))==0 %��Ϊ��б
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
    
else %��б
    if sum(ei(rmid-h,rmidc:n))==0 %��ת-44������������
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
            %����һ�Σ���h=1
            rmid=rmid+h; diff=ceil(diff/h); h=1;
        elseif rmidcol(length(rmidcol))>n
            if sum(ei(rmid,rmidcol(1):n))==0
                if rmid-1>0 && sum(ei(rmid-1,rmidcol(1):n))
                    rmid=rmid-1;
                    rmidc1=rmidcol(find(ei(rmid,rmidcol(1):n),1,'last'));
                else
                    %����һ�Σ���h=1
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
                %����һ�Σ���h=1
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

% %---------���½�---------
rmid=ceil(m/2); %�м���
rmidcol=ei(rmid,:);
h=10;
rmidc=find(rmidcol,1,'last');
rmid=rmid+h;
if sum(ei(rmid,rmidc+1:n))==0 %��Ϊ��б
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
    
else %��б
    if sum(ei(rmid+h,rmidc:n))==0 %��ת44������������
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
            %����һ�Σ���h=1
            rmid=rmid-h; diff=ceil(diff/h); h=1;
        elseif rmidcol(length(rmidcol))>n
            if sum(ei(rmid,rmidcol(1):n))==0
                if rmid+1<=m && sum(ei(rmid+1,rmidcol(1):n))
                    rmid=rmid+1;
                    rmidc1=rmidcol(find(ei(rmid,rmidcol(1):n),1,'last'));
                else
                    %����һ�Σ���h=1
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
                %����һ�Σ���h=1
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
