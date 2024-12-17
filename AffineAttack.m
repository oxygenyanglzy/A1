function img=AffineAttack(watermarkedimg,flag,b,c)
if flag==1
    % %------------仿射变换 -------%%%2019.12.09完成
    in_points = [1 b 0;c 1 0;0 0 1];%仿射变换
    %[0.9, 0.2, 0; 0.1, 1.2, 0; 0, 0, 1] 来源于对比文5：A robust image watermarking method based on DWT, DCT, and SVDusing a new technique for correction of main geometric attacksSaeid
    %     tform = maketform('affine', in_points); %创建结构体
    %     img = imtransform(watermarkedimg,tform);%执行仿射变换
    tform2=affine2d(in_points);
    img = imwarp(watermarkedimg,tform2);%执行仿射变换
%     figure(2), imshow(uint8(img)),title('仿射变换'); %显示变换结果
    %     imwrite(uint8(img),'Lena-Affine.jpg');
elseif flag==2
    %---------------------水平仿射变换--------------------%%
    in_points = [1 b 0;0 1 0;0 0 1];%水平仿射变换
    %     tform = maketform('affine',in_points); %创建新的恢复结构体2019-11-14 添加
    %     % tform = maketform('affine',[1 0 0; -0.5 1 0; 0 0 1]);
    %     img = imtransform(watermarkedimg,tform);
    tform2=affine2d(in_points);
    img = imwarp(watermarkedimg,tform2);%执行仿射变换
%     figure(2), imshow(uint8(img),'InitialMagnification','fit'),title('水平仿射变换');
    %     imwrite(uint8(img),'Lena-Horizontalshearing.jpg');
elseif flag==3
    %---------------------垂直仿射变换 --------------------%%
    in_points = [1 0 0;c 1 0;0 0 1];%垂直仿射变换
    tform = affine2d(in_points); %创建新的恢复结构体2019-11-14 添加
    % tform = maketform('affine',[1 0 0; -0.5 1 0; 0 0 1]);
    img = imwarp(watermarkedimg,tform);
    % tform2=projective2d(in_points);
    % img = imwarp(watermarkedimg,tform2);%执行仿射变换
%     figure(2), imshow(uint8(img)),title('垂直仿射变换');
    %     imwrite(uint8(img),'Lena-Verticalshearing.jpg');
end



