close all; clear; clc;
%%
tol         =   1e-12;
a           =   1;
%% func1
rho         =   1;
v           =   0;
tic;
I           =   HankelTransform(@(x)func1(x),rho,v,a,tol);
toc;
% correct answer    :   -1.000000000000000
fprintf('I1\t=\t%0.15f\n',I);
%% func2
rho         =   1;
v           =   1;
tic;
I           =   HankelTransform(@(x)func2(x),rho,v,a,tol);
toc;
% correct answer    :   +0.421024438240708
fprintf('I2\t=\t%0.15f\n',I);
%% func3
rho         =   1;
v           =   0;
tic;
I           =   HankelTransform(@(x)func3(x),rho,v,a,tol);
toc;
% correct answer    :   +1.000000000000000
fprintf('I3\t=\t%0.15f\n',I);
%% func4
rho         =   1;
v           =   10;
tic;
I           =   HankelTransform(@(x)func4(x),rho,v,a,tol);
toc;
% correct answer    :   +0.098970545308402
fprintf('I4\t=\t%0.15f\n',I);
%% func5
rho         =   200;
v           =   0;
tic;
I           =   HankelTransform(@(x)func5(x),rho,v,a,tol);
toc;
% correct answer    :   -1/rho^3
fprintf('I5\t=\t%0.15f\n',I);
%% Test Functions
function[y]=func1(x)
y               =   (x.^2);
end
function[y]=func2(x)
y               =   0.5*log(1+x.^2);
end
function[y]=func3(x)
y               =   (1-exp(-x))./(x.*log(1+sqrt(2)));
end
function[y]=func4(x)
y               =   x./(1+x.^2);
end
function[y]=func5(x)
y               =   x.^2;
end
%%
function[I]=HankelTransform(func,rho,v,a,tol)
%   function [I]=HankelTransform(func,rho,v,a,tol) computes
%   Hankel Transform for n order, based on Bessel Jn
%   In general HT = \int_{0}^{\inf} f(x) Jv(x\rho) dx
%
%  	func    :   function handle
%   rho     :   rho doesn't equal zeros
%   v       :   Bessel's order
%   a       :   first break point, defaul a=1
%   tol     :   tolerance
%%
switch nargin
    case 3
        a       =	1;
        tol     =	1e-12;
	case 4
        tol     =	1e-12;
end
b           =   FindBesselZero(a,rho,v);
q           =   pi/rho;
f           =   @(x)func(x).*besselj(v,rho*x);
I1        	=   quadgk(f,0,b);
I2        	=   PE(f,b,q,tol);
I          	=   I1+I2;
end
%%
function[y]=FindBesselZero(a,rho,v)
if rho<=0
    error('Incorrect input');
end
N_max           =   1000;
for i = 1:N_max
    B0          =   BesseljZeros(v,i);
    if B0 > a*rho
        y       =   B0/rho;
        break;
    end
end
end
%%
function[y]=BesseljZeros(v,n)
if n <= 10
    if v>10
        error('Bessel zeros not available. Please enter the first 10 zeros manually');
    end
    table(1,:)      =   [   2.404825557695773   5.520078110286311   8.653727912911013 ...
                            11.791534439014283  14.930917708487785  18.071063967910924 ...
                            21.211636629879258  24.352471530749302  27.493479132040253 ...
                            30.634606468431976  ];   
    table(2,:)      =   [   3.831705970207512   7.015586669815619  10.173468135062723 ...
                            13.323691936314223  16.470630050877634  19.615858510468243 ...
                            22.760084380592772  25.903672087618382  29.046828534916855 ...
                            32.189679910974405  ]; 
    table(3,:)      =   [   5.135622301840682   8.417244140399866  11.619841172149060 ...
                            14.795951782351260  17.959819494987826  21.116997053021844 ...
                            24.270112313573101  27.420573549984557  30.569204495516399 ...
                            33.716519509222699  ]; 
    table(4,:)      =   [   6.380161895923983   9.761023129981670  13.015200721698433 ...
                            16.223466160318768  19.409415226435012  22.582729593104443 ...
                            25.748166699294977  28.908350780921758  32.064852407097710 ...
                            35.218670738610115  ];
    table(5,:)      =   [   7.588342434503804  11.064709488501185  14.372536671617590 ...
                            17.615966049804832  20.826932956962388  24.019019524771110 ...
                            27.199087765981250  30.371007667117247  33.537137711819220 ...
                            36.699001128744648  ];
    table(6,:)      =   [   8.771483815959954  12.338604197466944  15.700174079711671 ...
                            18.980133875179920  22.217799896561267  25.430341154222702 ...
                            28.626618307291139  31.811716724047763  34.988781294559296 ...
                            38.159868561967130  ];
    table(7,:)      =   [   9.936109524217684  13.589290170541219  17.003819667816014 ...
                            20.320789213566506  23.586084435581391  26.820151983411403 ...
                            30.033722386570471  33.233041762847122  36.422019668258457 ...
                            39.603239416075404  ];
    table(8,:)      =   [   11.086370019245082  14.821268727013171  18.287582832481725 ...
                            21.641541019848400  24.934927887673023  28.191188459483200 ...
                            31.422794192265581  34.637089352069324  37.838717382853609 ...
                            41.030773691585537  ];
    table(9,:)      =   [   12.225092264004655  16.037774190887710  19.554536430997057 ...
                            22.945173131874618  26.266814641176641  29.545659670998546 ...
                            32.795800037341465  36.025615063869573  39.240447995178137 ...
                            42.443887743273557  ];
    table(10,:)      =   [  13.354300477435331  17.241220382489129  20.807047789264107 ...
                            24.233885257750551  27.583748963573004  30.885378967696674 ...
                            34.154377923855094  37.400099977156593  40.628553718964525 ...
                            43.843801420337350  ];
    table(11,:)      =   [  14.475500686554541  18.433463666966581  22.046985364697800 ...
                            25.509450554182823  28.887375063530460  32.211856199712727 ...
                            35.499909205373854  38.761807017881651  42.004190236671803 ...
                            45.231574103535046  ];
    y               =   table(v+1,n);
else 
    n               =   n-1;
    beta            =   (n+v/2-0.25)*pi;
    y               =   beta-(4*(v^2)-1)/(8*beta);
end
end
%%
function[y]=PE(func,a,q,tol)
kmax            =   25;
k               =   0;
while(k <= kmax)
    val = MosigMichalski(func,k,a,q);
    val(isnan(val))=0;
    val(isinf(abs(val)))=0;
    if k>1 && abs(val-old) < tol*abs(val) || abs(round(val,15))==0
        break;
    end
    old       	=   val;
    k         	=   k + 1;
    if k == kmax
        error('Maximum trials exceeded')
    end
end
y               =   val;
end
%%
function[y]=MosigMichalski(func,k,a,q)
mu                  =  	2;
R                   =   zeros(k+1,1);
eta                 =   zeros(k+1,1);
for n=0:k
    nn              =   n+1;
    R(nn,:)         =   Sn(func,n,a,q);
    eta(nn,:)       =   wn(func,n,a,q)./wn(func,n-1,a,q);
end
for n=1:k
    nn = n+1;
    for j=1:1:n
        d           =   xn(n-j+1,a,q) - xn(n-j,a,q);
        eta_k       =	eta(nn,:)./(1+mu*(j-1)*d/xn(n-j,a,q));
        R(nn-j,:)   =   (R(nn-j+1,:)-eta_k.*R(nn-j,:))./(1-eta_k);
    end
end
y = R(1,:);
end
%%
function [y] = Sn(func,n,a,q)
y               =	0;
for i=0:n
    y           =	y + u(func,i,a,q);
end
end
%%
function [y] = u(func,i,a,q)
y           =	quadgk(func,xn(i,a,q),xn(i+1,a,q));
end
%%
function [y] = wn(func,n,a,q)
y           =	u(func,n,a,q);  
end
%%
function [y] = xn(n,a,q)
y               =	a + n*q;
end
%%
