function afterAffine= Affine(Image,~,~,~,~,Frequency,crypt)
                         %ͼ����ֵ�������任ͼ�����Һ���
%�������
%		Image:			   �����ܣ�������ܣ�ͼ���ļ�����ע��д��ʽ��׺������СΪ��ά
%       a,b,c,d:             ��Կ       ����任�ľ���
%		Frequency:		 ͼ����Ҫ�任�����Ĵ���
%       crypt:               0������;    1������
%�������
%	afterAffine:	     ת�����ͼ�����ݾ���
%                             �����M��Ӧ��ͼ���ļ�
%  2019.1.13  ���

SizeImg=size(Image);
%���ͼ���Ƕ�ά���򲻴�������
if (length(SizeImg) == 2)
    if SizeImg(1) ~= SizeImg(2)
        disp('���Ƿ��󣬲��ܷ���任��');
        return
    end
else
    disp('���Ƕ�ά���飬���ܽ��з���任��');
    return
end

N=SizeImg(1);   %NΪ����ĳ��ȣ����ȣ�
P=[1,-1;-1,0];
if crypt==0       %���ܹ���
    img1=Image;
    for t=1:Frequency
        for x=1:N        %��Ϊ�����±�û��0�±꣬�����Ǵ�1��N�������Ǵ�0��N-1��
            for y=1:N
                if x<=y-1         %�����ж�����Ϊx<=y-1��������x<=y��
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
elseif crypt==1        %���ܹ���
    Q=inv(P);             %Q����ΪP����������
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
    disp('���ܻ���ܲ�������!  ��ʾ��0��ʾ���ܣ�1��ʾ���ܣ�');
end
clear P;
end

