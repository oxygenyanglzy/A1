function img=AffineAttack(watermarkedimg,flag,b,c)
if flag==1
    % %------------����任 -------%%%2019.12.09���
    in_points = [1 b 0;c 1 0;0 0 1];%����任
    %[0.9, 0.2, 0; 0.1, 1.2, 0; 0, 0, 1] ��Դ�ڶԱ���5��A robust image watermarking method based on DWT, DCT, and SVDusing a new technique for correction of main geometric attacksSaeid
    %     tform = maketform('affine', in_points); %�����ṹ��
    %     img = imtransform(watermarkedimg,tform);%ִ�з���任
    tform2=affine2d(in_points);
    img = imwarp(watermarkedimg,tform2);%ִ�з���任
%     figure(2), imshow(uint8(img)),title('����任'); %��ʾ�任���
    %     imwrite(uint8(img),'Lena-Affine.jpg');
elseif flag==2
    %---------------------ˮƽ����任--------------------%%
    in_points = [1 b 0;0 1 0;0 0 1];%ˮƽ����任
    %     tform = maketform('affine',in_points); %�����µĻָ��ṹ��2019-11-14 ���
    %     % tform = maketform('affine',[1 0 0; -0.5 1 0; 0 0 1]);
    %     img = imtransform(watermarkedimg,tform);
    tform2=affine2d(in_points);
    img = imwarp(watermarkedimg,tform2);%ִ�з���任
%     figure(2), imshow(uint8(img),'InitialMagnification','fit'),title('ˮƽ����任');
    %     imwrite(uint8(img),'Lena-Horizontalshearing.jpg');
elseif flag==3
    %---------------------��ֱ����任 --------------------%%
    in_points = [1 0 0;c 1 0;0 0 1];%��ֱ����任
    tform = affine2d(in_points); %�����µĻָ��ṹ��2019-11-14 ���
    % tform = maketform('affine',[1 0 0; -0.5 1 0; 0 0 1]);
    img = imwarp(watermarkedimg,tform);
    % tform2=projective2d(in_points);
    % img = imwarp(watermarkedimg,tform2);%ִ�з���任
%     figure(2), imshow(uint8(img)),title('��ֱ����任');
    %     imwrite(uint8(img),'Lena-Verticalshearing.jpg');
end



