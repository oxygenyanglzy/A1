% Test code to verify DHT and IDHT
A = [196,20,33,175;
    37,56,79,188;
    100,27,89,33;
    27,55,36,99];
% Define the radial function for dht
h = @(r) A;
% Parameters.
R =3;  % Maximum radius
N = 4;  % Number of sampling points
n = 0;  % Transform order
% Perform the DHT
[H1,k,r,I1,K1,R1,h]=dht(h,R,N,n);
H2=[948.75,487.236550241111,785.669936634047,1377.93590326907;22.9388804601482,-160.592089312480,-194.655068041473,34.3175214367763;114.553085253851,63.2790352094624,-81.2928816692196,159.345975397545;209.989466293873,-103.113030389414,23.0416055717673,-277.264286518022];
recover = idht(H2, I1, K1, R1);


recover = int16(recover);

