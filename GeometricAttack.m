function [attackimg,aname]=GeometricAttack(watermarkedimg)
para=20;attackimg=RotatingAttack(watermarkedimg,para); name='rotate';%��ת
% para=1;tar1=[50,1;512,60;450,512;1,400]; attackimg=SlantAttack(watermarkedimg,tar1); name='slant'; %��б
% para=2;tar2=[1,120;350,1;512,420;100,512]; attackimg=SlantAttack(watermarkedimg,tar2); name='slant';%��б
% para=0.09;attackimg=AffineAttack(watermarkedimg,1,para,para); name='affine';% ����
% para=0.1;attackimg=AffineAttack(watermarkedimg,2,para,0); name='y-affine';% ˮƽ���� ��Ե�����Ҳ�����
% para=0.1;attackimg=AffineAttack(watermarkedimg,3,0,para); name='x-affine';% ��ֱ����
% para=-1;attackimg=TranslatingAttack(watermarkedimg,-35,-25); name='translate';% ƽ�� ��Ե�����Ҳ�����
% imshow(attackimg);
aname=[name,num2str(para)]; %rotate30
