function M=Arnold(Image,Frequency,crypt)
%图像数值矩阵Arnold转换函数
%输入参数
%		Image:			待加密（待解密）图像文件名（注意写格式后缀），只能为二维
%		Frequency:		图像需要变换迭的次数
%       crypt           0～加密;1～解密

%输出参数
%		M:				转换后图像数据矩阵
%                       输出由M对应的图像文件
if nargin<3
    disp('请按程序的输入参数格式输入参数！！！');
    return;
end

if crypt~=0 & crypt~=1
    disp('encrypt 必须为0或1！');
end

%将Q赋值给M，计算Q的大小
%Q=imread(Image);
Q=Image;
M = Q ;
Size_Q   = size(Q);

%如果不是二维或三维数组，则不处理，返回
if (length(Size_Q) == 2) 
   if Size_Q(1) ~= Size_Q(2)  
      disp('不是方阵，不能Arnold转换');
      return
   end
else
   disp('不是二维数组，不进行Arnold变换');
   return 
end

    %------------------------------------------
   %Arnold转换
   n = 0;
   K = Size_Q(1);
   
   M1_t = Q;
   M2_t = Q;
   
   if crypt==1   %解密
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
% 求周期
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