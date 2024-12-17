function M=Arnold(Image,Frequency,crypt)
%ͼ����ֵ����Arnoldת������
%�������
%		Image:			�����ܣ������ܣ�ͼ���ļ�����ע��д��ʽ��׺����ֻ��Ϊ��ά
%		Frequency:		ͼ����Ҫ�任���Ĵ���
%       crypt           0������;1������

%�������
%		M:				ת����ͼ�����ݾ���
%                       �����M��Ӧ��ͼ���ļ�
if nargin<3
    disp('�밴��������������ʽ�������������');
    return;
end

if crypt~=0 & crypt~=1
    disp('encrypt ����Ϊ0��1��');
end

%��Q��ֵ��M������Q�Ĵ�С
%Q=imread(Image);
Q=Image;
M = Q ;
Size_Q   = size(Q);

%������Ƕ�ά����ά���飬�򲻴�������
if (length(Size_Q) == 2) 
   if Size_Q(1) ~= Size_Q(2)  
      disp('���Ƿ��󣬲���Arnoldת��');
      return
   end
else
   disp('���Ƕ�ά���飬������Arnold�任');
   return 
end

    %------------------------------------------
   %Arnoldת��
   n = 0;
   K = Size_Q(1);
   
   M1_t = Q;
   M2_t = Q;
   
   if crypt==1   %����
       Frequency=ArnoldPeriod( Size_Q(1) )-Frequency;
   end
       
   for s = 1:Frequency
       n = n + 1;
       if mod(n,2) == 0
            for i = 1:K
               for j = 1:K
                  c = M2_t(i,j);
                  M1_t(mod(i+j-2,K)+1,mod(i+2*j-3,K)+1) = c;
               end
            end
       else
            for i = 1:K
               for j = 1:K
                   c = M1_t(i,j);
                   M2_t(mod(i+j-2,K)+1,mod(i+2*j-3,K)+1) = c;
               end
            end
       end
   end
   
   if mod(Frequency,2) == 0
      M = M1_t;
   else
      M = M2_t;
   end
   %------------------------------------------
   %imwrite( double(M)/255,strcat( 'Arnold_',num2str(Frequency),'_',Image ),'bmp' );
   %imshow(M);
   
function Period=ArnoldPeriod(N)
% ������
if ( N<2 )
    Period=0;
    return;
end

n=1;
x=1;
y=1;
while n~=0
    xn=x+y;
    yn=x+2*y;
    if ( mod(xn,N)==1 && mod(yn,N)==1 )
        Period=n;
        return;
    end
    x=mod(xn,N);
    y=mod(yn,N);
    n=n+1;
end