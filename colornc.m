%Name:		wei--kai19830426@163.com
%Course:	Digital Watermarking
%Project: 	Calculates the NC (Normalized Correlation)
%function:  of watermark images A and A', both of size MxN
%%%%%%% NC %%%%%%

function N=colornc(mark_get,mark_prime)
nclevel=0;
mark_get=double(mark_get);
mark_prime=double(mark_prime);
if size(mark_get)~=size(mark_prime)
    error('Input vectors must  be the same size!')
else
    [m,n]=size(mark_get);
    fenzi=0;
    fenmu1=0;
    fenmu2=0;
    for level=1:3
    for i=1:m
        for j=1:m
            fenzi=fenzi+mark_get(i,j,level)*mark_prime(i,j,level);
            fenmu1=fenmu1+mark_prime(i,j,level)*mark_prime(i,j,level);
            fenmu2=fenmu2+mark_get(i,j,level)*mark_get(i,j,level);
        end
    end
    ss=min(fenzi/fenmu1,fenzi/fenmu2);
    nclevel=nclevel+ss;
    end
    N=nclevel/3;
end