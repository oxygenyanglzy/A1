function img=TranslatingAttack(watermarkedimg,x,y)
%----------ƽ�Ʊ任---------------2019.12.06���
I=watermarkedimg;
se = translate(strel(1), [x y]);
img = imdilate(I,se);
% imwrite(uint8(img),['Lena-Translation(',num2str(x),',',num2str(y),').jpg']);
% % imwrite(uint8(I1),'ImageDatebases\Lena-Translation(-20,0).jpg');
% figure(2),imshow(uint8(img),'InitialMagnification','fit');%ͼ�����ʱ��ʹ�ø�����Ӧ������ʾ
end