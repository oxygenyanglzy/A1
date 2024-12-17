function afterAffine= Affine(Image,~,~,~,~,Frequency,crypt)
                         %图像数值矩阵仿射变换图像置乱函数
%输入参数
%		Image:			   待加密（或待解密）图像文件名（注意写格式后缀），大小为二维
%       a,b,c,d:             秘钥       仿射变换的矩阵
%		Frequency:		 图像需要变换更迭的次数
%       crypt:               0～加密;    1～解密
%输出参数
%	afterAffine:	     转换后的图像数据矩阵
%                             输出由M对应的图像文件
%  2019.1.13  完成

SizeImg=size(Image);
%如果图像不是二维，则不处理，返回
if (length(SizeImg) == 2)
    if SizeImg(1) ~= SizeImg(2)
        disp('不是方阵，不能仿射变换！');
        return
    end
else
    disp('不是二维数组，不能进行仿射变换！');
    return
end

N=SizeImg(1);   %N为方阵的长度（或宽度）
P=[1,-1;-1,0];
if crypt==0       %加密过程
    img1=Image;
    for t=1:Frequency
        for x=1:N        %因为矩阵下标没有0下标，所以是从1到N，而不是从0到N-1；
            for y=1:N
                if x<=y-1         %进而判断条件为x<=y-1，而不是x<=y；
                    xx=P(1,1)*x+P(1,2)*y+N+1;
                    yy=P(2,1)*x+P(2,2)*y+N+1;
                else
                    xx=P(1,1)*x+P(1,2)*y+1;
                    yy=P(2,1)*x+P(2,2)*y+N+1;
                end
                img2(xx,yy)=img1(x,y);
            end
        end
        img1=img2;
    end
    afterAffine=img2;
elseif crypt==1        %解密过程
    Q=inv(P);             %Q矩阵为P矩阵的逆矩阵
    img3=Image;
    for t=1:Frequency
        for x=1:N
            for y=1:N
                if x+y<=N+1
                    xx=Q(1,1)*x+Q(1,2)*y+N+1;
                    yy=Q(2,1)*x+Q(2,2)*y+N+2;
                else
                    xx=Q(1,1)*x+Q(1,2)*y+N+1;
                    yy=Q(2,1)*x+Q(2,2)*y+2*N+2;
                end
                img4(xx,yy)=img3(x,y);
            end
        end
        img3=img4;
    end
    afterAffine=img4;
else
    disp('加密或解密参数错误!  提示：0表示加密，1表示解密！');
end
clear P;
end

