A=[143,164,152,133;129,167,177,158;129,153,172,177;131,146,158,158];
[L,U]=lu(A);
A1=[128,115,113,89;28,56,90,1;25,45,25,55;184,32,0,15];
[L1,U1]=lu(A1);
H='lena.jpg';
Host=imread(H);%��������ͼƬ
hfx=Host(:,:,1);
% I = rgb2gray(Host);
A2=[225	224	224	223;
224	224	224	224;
224	224	225	223;
224	225	225	223];
[L2,U2]=lu(A2);
h=sum(U2(1,:));
a1=[49,46,0,0;43,42,0,46;41,42,55,50;42,44,57,50];
[l1,u1]=lu(a1)
% a2=l1*u1;
% u2=u1;
% for n=1:4
%     u2(1,n)=u2(1,n)+10;
% end
% a2=l1*u2;
% [l2,u3]=lu(a2);
