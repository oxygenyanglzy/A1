close all;
clear;
clc;
Hdata=["pine.jpg"];
Wdata=["qq.jpg"];
count=1;
for Hi=1:1
    H=Hdata(Hi);
    host =imread(H);
    for Wi=1:1
        ending=2;
        %athens.jpg；house.jpg；sailboat.jpg；peppers.jpg；； butrfly1.jpg；sea.jpg；bluheron.jpg；couple.jpg；
        W=Wdata(Wi);   % ldu.jpg；200100732.jpg；
        % T=100;   %T=10
        watermarkimage=imread(W);
        waterlenth=size(watermarkimage,1)*size(watermarkimage,2)*8;
        for level=1:3
            watermarklevel=watermarkimage(:,:,level);
            bitwater(level,:)=gainlevelwatermark(watermarklevel);
            H=double(host(:,:,level));
            watermark=reshape(bitwater(level,:),128,64);
            %                 [nrw,ncw,vim]=size(watermark);
            %                 [nrf,ncf]=size(host(:,:,level));
            %                 rng(nrw);
            %                 pn_key_r=randperm(floor(nrf/4),nrw);
            %                 rand_r=pn_key_r;%%选取嵌入水印的位置
            %                 rng(ncw);
            %                 pn_key_c=randperm(floor(ncf/4),ncw);
            %                 rand_c=pn_key_c;
            %参数的设置
            lambda=216; D=0.032;
            delta=0.5*D;
            psi=D;  %这个参数是啥？？
            rho=1.5*D; eta=2*D;
            blockSize=4;
            for k=1:128%%水印图像中0和1的个数
                for l=1:64
                    A=H((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
                    A_average=mean(A(:));
                    A_lambda=A+(lambda-A_average)*ones(blockSize)+(rand(blockSize)-0.5);%Level shifting with dither noise
                    [u,s,v]=svd(A_lambda);
                    rc=u(:,1)'*A_lambda*ones(blockSize,1)+ones(1,blockSize)*A_lambda*v(:,1);  %sign correction
                    if rc<0
                        u(:,1)=-u(:,1);
                        v(:,1)=-v(:,1);
                    end
                    %mixed modulation
                    g=u(2,1)-u(3,1);
                    if watermark(k,l)=='1'
                        gamma1=min(psi, max(g,delta));
                    else
                        gamma1=max(-psi, min(g,-delta));
                    end
                    if g>rho
                        if watermark(k,l)=='1'
                            gamma21=max(0, floor((g-eta)/(2*D))*2*D+D+eta);
                        else
                            gamma21=floor((g-eta)/(2*D)+0.5)*2*D+eta;
                        end
                    elseif g<-rho
                        if watermark(k,l)=='1'
                            gamma22=ceil((g+eta)/(2*D)-0.5)*2*D-eta;
                        else
                            gamma22=min(0, ceil((g+eta)/(2*D))*2*D-D-eta);
                        end
                    end
                    if g>rho && abs(gamma21-g)<abs(gamma1-g)
                        g_hat=gamma21;
                        zta1=0.5*(-u(2,1)+g+u(3,1));  % 文中错误的计算公式 zta1=0.5*g_hat-u(2,1);
                        zta2=0.5*(g_hat+u(3,1)-u(2,1));
                        omega=0.5; zta=omega*zta1+(1-omega)*zta2; %Orthonormal restoration  zta=1*zta1;
                        u_hat=u;     u_hat(2,1)=u(2,1)+zta;           u_hat(3,1)=u(2,1)-g_hat+zta;
                        for kk = [1, blockSize]
                            u_hat(kk,1)=u(kk,1)*sqrt((1-u_hat(2,1)^2-u_hat(3,1)^2)/(1-u(2,1)^2-u(3,1)^2));
                        end
                    elseif g<-rho && abs(gamma22-g)<abs(gamma1-g)
                        g_hat=gamma22;
                        zta1=0.5*(-u(2,1)+g+u(3,1));
                        zta2=0.5*(g_hat+u(3,1)-u(2,1));
                        omega=0.5;  zta=omega*zta1+(1-omega)*zta2; %   zta=1*zta1;
                        u_hat=u;     u_hat(2,1)=u(2,1)+zta;           u_hat(3,1)=u(2,1)-g_hat+zta;
                        for kk = [1, blockSize]
                            u_hat(kk,1)=u(kk,1)*sqrt((1-u_hat(2,1)^2-u_hat(3,1)^2)/(1-u(2,1)^2-u(3,1)^2));
                        end
                    else
                        g_hat=gamma1;
                        zta1=0.5*(-u(2,1)+g+u(3,1));
                        zta2=0.5*(g_hat+u(3,1)-u(2,1));
                        omega=0.25; zta=omega*zta1+(1-omega)*zta2; %zta=1*zta1;
                        u_hat=u;     u_hat(2,1)=u(2,1)+zta;           u_hat(3,1)=u(2,1)-g_hat+zta;
                        for kk = [2,3]
                            u_hat(kk,1)=u_hat(kk,1)*sqrt((1-u(2,1)^2-u(3,1)^2)/(1-u_hat(2,1)^2-u_hat(3,1)^2));
                        end
                    end
                    u_hat=gschmidt(u_hat); %%斯密特正交化
                    s_hat=s;%失真补偿
                    for kk = 1:blockSize
                        P=u_hat(:,kk)*v(:,kk)';
                        temp=0;
                        for i = 1: blockSize
                            temp=temp+A_lambda(i,:)*P(i,:)';
                        end
                        s_hat(kk,kk)=temp;
                    end
                    A_hat=u_hat*s_hat*v'+(A_average-lambda)*ones(blockSize,blockSize);
                    A_hat=uint8(round(A_hat));
                    H((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize)=A_hat;
                end
            end
            rgb_image(:,:,level)=H;
        end
        imwrite(uint8(rgb_image),'image/WatermarkedImage/watermarked.bmp');
        psnrval=colorpsnr(host,rgb_image);
        ssimval=0;
        for level=1:3
            ssimval=ssimval+colorssim(host(:,:,level), rgb_image(:,:,level));
        end
        ssimval=ssimval/3;
        %             figure();
        %             subplot(121),imshow('WatermarkedImage/watermarked.bmp'),title([' PSNR/SSIM=',num2str(psnrval),'/',num2str(ssimval)]);
        
        watermarkedimage = imread('image/WatermarkedImage/watermarked.bmp');
        for level=1:3
            watermarklevel=watermarkimage(:,:,level);
            bitwater(level,:)=gainlevelwatermark(watermarklevel);
            Hd=double(watermarkedimage(:,:,level));
            watermark=reshape(bitwater(level,:),128,64);
            %                 [nrw,ncw,vim]=size(watermark);
            %                 [nrf,ncf]=size(host(:,:,level));
            %                 rng(nrw);
            %                 pn_key_r=randperm(floor(nrf/4),nrw);
            %                 rand_r=pn_key_r;%%选取嵌入水印的位置
            %                 rng(ncw);
            %                 pn_key_c=randperm(floor(ncf/4),ncw);
            %                 rand_c=pn_key_c;
            %参数的设置
            %                 lambda=216; D=0.032;
            %                 delta=0.5*D;
            %                 psi=D;  %这个参数是啥？？
            %                 rho=1.5*D; eta=2*D;
            %                 blockSize=4;
            iter=1;
            lambda=216; D=0.032;
            delta=0.5*D;
            psi=D;
            rho=1.5*D; eta=2*D;
            blockSize=4;
            % NumOfColumn=256;
            % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
            for k=1:128
                for l=1:64
                    DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
                    DA_average=mean(DA(:));
                    DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
                    [ud,sd,vd]=svd(DA_lambda);
                    rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
                    if rc<0
                        ud(:,1)=-ud(:,1);
                        vd(:,1)=-vd(:,1);
                    end
                    g_hat=ud(2,1)-ud(3,1);
                    if abs(g_hat)<=rho
                        if g_hat>=0
                            Dw(k,l)='1';
                        else
                            Dw(k,l)='0';
                        end
                    else
                        if g_hat>rho
                            Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
                        elseif g_hat<-rho
                            Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
                        end
                    end
                end
            end
            wa=watermark(:)';
            p=0;
            for i=1:128
                for j=1:64
                    if Dw(i,j)~=watermark(i,j)
                        p=p+1;
                    end
                end
            end
            ber(level)=p/(8192);
            Dw2=reshape(Dw,1,8192);
            dd='';
            for i=1:8192
                dd=strcat(dd,num2str(Dw2(i)));
            end
            extractwatermarkimage(:,:,level)=gainlevelimage(dd);
            %                 watermark_d=reshape(Dw,[128,64]);
        end
        imwrite(uint8(extractwatermarkimage),'image/ext/extractwater.bmp');
        ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
        berval=mean(ber);
        result(count,ending)=psnrval;
        ending=ending+1;
        result(count,ending)=ssimval;
        ending=ending+1;
        result(count,ending)=ncval;
        ending=ending+1;
        result(count,ending)=berval;
        %             time1(count,ending)=time1;
        %             time2(count,ending)=time2;
        ending=ending+1;
        figure();
        subplot(121),imshow('WatermarkedImage/watermarked.bmp'),title([' PSNR/SSIM=',num2str(psnrval),'/',num2str(ssimval)]);
        subplot(122),imshow('image/ext/extractwater.bmp'),title([' NC/BER=',num2str(ncval),'/',num2str(berval),]);
        % % JPEG测试，分层攻击图像
%         watermarkedimage = imread('image/WatermarkedImage/watermarked.bmp');
%         for Q=50:10:90
%             for level=1:3
%                 imwrite(watermarkedimage(:,:,level),'image/WatermarkedImage/JPEG.bmp','jpg','Quality',Q);
%                 watermarklevel=watermarkimage(:,:,level);
%                 bitwater(level,:)=gainlevelwatermark(watermarklevel);
%                 Hd=double(imread('image/WatermarkedImage/JPEG.bmp'));   %temp.bmp 为压缩后的图片
%                 watermark=reshape(bitwater(level,:),128,64);
%                 %                 [nrw,ncw,vim]=size(watermark);
%                 %                 [nrf,ncf]=size(host(:,:,level));
%                 %                 rng(nrw);
%                 %                 pn_key_r=randperm(floor(nrf/4),nrw);
%                 %                 rand_r=pn_key_r;%%选取嵌入水印的位置
%                 %                 rng(ncw);
%                 %                 pn_key_c=randperm(floor(ncf/4),ncw);
%                 %                 rand_c=pn_key_c;
%                 %参数的设置
%                 %                 lambda=216; D=0.032;
%                 %                 delta=0.5*D;
%                 %                 psi=D;  %这个参数是啥？？
%                 %                 rho=1.5*D; eta=2*D;
%                 %                 blockSize=4;
%                 iter=1;
%                 lambda=216; D=0.032;
%                 delta=0.5*D;
%                 psi=D;
%                 rho=1.5*D; eta=2*D;
%                 blockSize=4;
%                 % NumOfColumn=256;
%                 % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
%                 for k=1:128
%                     for l=1:64
%                         DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
%                         DA_average=mean(DA(:));
%                         DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
%                         [ud,sd,vd]=svd(DA_lambda);
%                         rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
%                         if rc<0
%                             ud(:,1)=-ud(:,1);
%                             vd(:,1)=-vd(:,1);
%                         end
%                         g_hat=ud(2,1)-ud(3,1);
%                         if abs(g_hat)<=rho
%                             if g_hat>=0
%                                 Dw(k,l)='1';
%                             else
%                                 Dw(k,l)='0';
%                             end
%                         else
%                             if g_hat>rho
%                                 Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
%                             elseif g_hat<-rho
%                                 Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
%                             end
%                         end
%                     end
%                 end
%                 wa=watermark(:)';
%                 p=0;
%                 for i=1:128
%                     for j=1:64
%                         if Dw(i,j)~=watermark(i,j)
%                             p=p+1;
%                         end
%                     end
%                 end
%                 ber(level)=p/(8192);
%                 Dw2=reshape(Dw,1,8192);
%                 dd='';
%                 for i=1:8192
%                     dd=strcat(dd,num2str(Dw2(i)));
%                 end
%                 extractwatermarkimage(:,:,level)=gainlevelimage(dd);
%                 %                 watermark_d=reshape(Dw,[128,64]);
% 
%             end
%             % for i=1:32
%             %  error_rate = calculate_error_rate(extractwatermarkimage(i,:,3),watermarkimage(i,:,3));
%             %   disp(['Error rate: ' num2str(error_rate)]);
%             % end
%             imwrite(uint8(extractwatermarkimage),'image/ext/extractwater_JPEG.bmp');
%             ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
%             berval=mean(ber);
%             result(count,ending)=ncval;
%             ending=ending+1;
%             result(count,ending)=berval;
%             %             time1(count,ending)=time1;
%             %             time2(count,ending)=time2;
%             ending=ending+1;
% %                             figure();
% %                             subplot(121),imshow('WatermarkedImage/watermarked.bmp'),title([' PSNR/SSIM=',num2str(psnrval),'/',num2str(ssimval)]);
% %                             subplot(122),imshow('image/ext/extractwater_JPEG.bmp'),title([' NC/BER=',num2str(ncval),'/',num2str(berval),]);
%             %         figure();
%             %         subplot(121),imshow('WatermarkedImage/watermarked.bmp'),title([num2str(Q),'   PSNR=',num2str(psnrval)]);
%             %         subplot(122),imshow('image/ext/extractwater_JPEG.bmp'),title([' NC/BER=',num2str(ncval),'/',num2str(berval),]);
%         end
%         % % JPEG2000测试，分层攻击图像
%         watermarkedimage = imread("image/WatermarkedImage/watermarked.bmp");
%         for  Q=4:1:6
%             for level=1:3
%                 imwrite(watermarkedimage(:,:,level),'image/WatermarkedImage/JPEG2000.bmp','jp2','compressionratio',Q);
%                 hostlevel=double(imread('image/WatermarkedImage/JPEG2000.bmp'));   %temp.bmp 为压缩后的图片
%                 watermarklevel=watermarkimage(:,:,level);
%                 bitwater(level,:)=gainlevelwatermark(watermarklevel);
%                 Hd=hostlevel;   %temp.bmp 为压缩后的图片
%                 watermark=reshape(bitwater(level,:),128,64);
%                 %                 [nrw,ncw,vim]=size(watermark);
%                 %                 [nrf,ncf]=size(host(:,:,level));
%                 %                 rng(nrw);
%                 %                 pn_key_r=randperm(floor(nrf/4),nrw);
%                 %                 rand_r=pn_key_r;%%选取嵌入水印的位置
%                 %                 rng(ncw);
%                 %                 pn_key_c=randperm(floor(ncf/4),ncw);
%                 %                 rand_c=pn_key_c;
%                 %参数的设置
%                 %                 lambda=216; D=0.032;
%                 %                 delta=0.5*D;
%                 %                 psi=D;  %这个参数是啥？？
%                 %                 rho=1.5*D; eta=2*D;
%                 %                 blockSize=4;
%                 iter=1;
%                 lambda=216; D=0.032;
%                 delta=0.5*D;
%                 psi=D;
%                 rho=1.5*D; eta=2*D;
%                 blockSize=4;
%                 % NumOfColumn=256;
%                 % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
%                 for k=1:128
%                     for l=1:64
%                         DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
%                         DA_average=mean(DA(:));
%                         DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
%                         [ud,sd,vd]=svd(DA_lambda);
%                         rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
%                         if rc<0
%                             ud(:,1)=-ud(:,1);
%                             vd(:,1)=-vd(:,1);
%                         end
%                         g_hat=ud(2,1)-ud(3,1);
%                         if abs(g_hat)<=rho
%                             if g_hat>=0
%                                 Dw(k,l)='1';
%                             else
%                                 Dw(k,l)='0';
%                             end
%                         else
%                             if g_hat>rho
%                                 Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
%                             elseif g_hat<-rho
%                                 Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
%                             end
%                         end
%                     end
%                 end
%                 wa=watermark(:)';
%                 p=0;
%                 for i=1:128
%                     for j=1:64
%                         if Dw(i,j)~=watermark(i,j)
%                             p=p+1;
%                         end
%                     end
%                 end
%                 ber(level)=p/(8192);
%                 Dw2=reshape(Dw,1,8192);
%                 dd='';
%                 for i=1:8192
%                     dd=strcat(dd,num2str(Dw2(i)));
%                 end
%                 extractwatermarkimage(:,:,level)=gainlevelimage(dd);
%                 %                 watermark_d=reshape(Dw,[128,64]);
% 
%             end
%             % for i=1:32
%             %  error_rate = calculate_error_rate(extractwatermarkimage(i,:,3),watermarkimage(i,:,3));
%             %   disp(['Error rate: ' num2str(error_rate)]);
%             % end
%             imwrite(uint8(extractwatermarkimage),'image/ext/extractwater_JPEG.bmp');
%             ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
%             berval=mean(ber);
%             result(count,ending)=ncval;
%             ending=ending+1;
%             result(count,ending)=berval;
%             %             time1(count,ending)=time1;
%             %             time2(count,ending)=time2;
%             ending=ending+1;
%         end
%         % % 中值滤波攻击测试，分层攻击图像
        watermarkedimage = imread("image/WatermarkedImage/watermarked.bmp");
        for Q=3:2:5
            for level=1:3
                filtering=medfilt2(watermarkedimage(:,:,level),[Q,Q]); %二维中值滤波
                watermarklevel=watermarkimage(:,:,level);
                bitwater(level,:)=gainlevelwatermark(watermarklevel);
                Hd=double(filtering);   %temp.bmp 为压缩后的图片
                watermark=reshape(bitwater(level,:),128,64);
                iter=1;
                lambda=216; D=0.032;
                delta=0.5*D;
                psi=D;
                rho=1.5*D; eta=2*D;
                blockSize=4;
                % NumOfColumn=256;
                % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
                for k=1:128
                    for l=1:64
                        DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
                        DA_average=mean(DA(:));
                        DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
                        [ud,sd,vd]=svd(DA_lambda);
                        rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
                        if rc<0
                            ud(:,1)=-ud(:,1);
                            vd(:,1)=-vd(:,1);
                        end
                        g_hat=ud(2,1)-ud(3,1);
                        if abs(g_hat)<=rho
                            if g_hat>=0
                                Dw(k,l)='1';
                            else
                                Dw(k,l)='0';
                            end
                        else
                            if g_hat>rho
                                Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
                            elseif g_hat<-rho
                                Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
                            end
                        end
                    end
                end
                wa=watermark(:)';
                p=0;
                for i=1:128
                    for j=1:64
                        if Dw(i,j)~=watermark(i,j)
                            p=p+1;
                        end
                    end
                end
                ber(level)=p/(8192);
                Dw2=reshape(Dw,1,8192);
                dd='';
                for i=1:8192
                    dd=strcat(dd,num2str(Dw2(i)));
                end
                extractwatermarkimage(:,:,level)=gainlevelimage(dd);
                %                 watermark_d=reshape(Dw,[128,64]);

            end
            % for i=1:32
            %  error_rate = calculate_error_rate(extractwatermarkimage(i,:,3),watermarkimage(i,:,3));
            %   disp(['Error rate: ' num2str(error_rate)]);
            % end
            imwrite(uint8(extractwatermarkimage),'image/ext/extractwater_JPEG.bmp');
            ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
            berval=mean(ber);
            result(count,ending)=ncval;
            ending=ending+1;
            result(count,ending)=berval;
            %             time1(count,ending)=time1;
            %             time2(count,ending)=time2;
            ending=ending+1;
        end
%         % % 高斯噪声攻击测试，分层攻击图像
%         watermarkedimage = imread("image/WatermarkedImage/watermarked.bmp");
%         for Q=3:1:5
%             for level=1:3
%                 filtering=imfilter(watermarkedimage(:,:,level),fspecial('gaussian',[Q,Q])); %二维中值滤波
%                 watermarklevel=watermarkimage(:,:,level);
%                 bitwater(level,:)=gainlevelwatermark(watermarklevel);
%                 Hd=double(filtering);   %temp.bmp 为压缩后的图片
%                 watermark=reshape(bitwater(level,:),128,64);
%                 %                 [nrw,ncw,vim]=size(watermark);
%                 %                 [nrf,ncf]=size(host(:,:,level));
%                 %                 rng(nrw);
%                 %                 pn_key_r=randperm(floor(nrf/4),nrw);
%                 %                 rand_r=pn_key_r;%%选取嵌入水印的位置
%                 %                 rng(ncw);
%                 %                 pn_key_c=randperm(floor(ncf/4),ncw);
%                 %                 rand_c=pn_key_c;
%                 %参数的设置
%                 %                 lambda=216; D=0.032;
%                 %                 delta=0.5*D;
%                 %                 psi=D;  %这个参数是啥？？
%                 %                 rho=1.5*D; eta=2*D;
%                 %                 blockSize=4;
%                 iter=1;
%                 lambda=216; D=0.032;
%                 delta=0.5*D;
%                 psi=D;
%                 rho=1.5*D; eta=2*D;
%                 blockSize=4;
%                 % NumOfColumn=256;
%                 % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
%                 for k=1:128
%                     for l=1:64
%                         DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
%                         DA_average=mean(DA(:));
%                         DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
%                         [ud,sd,vd]=svd(DA_lambda);
%                         rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
%                         if rc<0
%                             ud(:,1)=-ud(:,1);
%                             vd(:,1)=-vd(:,1);
%                         end
%                         g_hat=ud(2,1)-ud(3,1);
%                         if abs(g_hat)<=rho
%                             if g_hat>=0
%                                 Dw(k,l)='1';
%                             else
%                                 Dw(k,l)='0';
%                             end
%                         else
%                             if g_hat>rho
%                                 Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
%                             elseif g_hat<-rho
%                                 Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
%                             end
%                         end
%                     end
%                 end
%                 wa=watermark(:)';
%                 p=0;
%                 for i=1:128
%                     for j=1:64
%                         if Dw(i,j)~=watermark(i,j)
%                             p=p+1;
%                         end
%                     end
%                 end
%                 ber(level)=p/(8192);
%                 Dw2=reshape(Dw,1,8192);
%                 dd='';
%                 for i=1:8192
%                     dd=strcat(dd,num2str(Dw2(i)));
%                 end
%                 extractwatermarkimage(:,:,level)=gainlevelimage(dd);
%                 %                 watermark_d=reshape(Dw,[128,64]);
% 
%             end
%             % for i=1:32
%             %  error_rate = calculate_error_rate(extractwatermarkimage(i,:,3),watermarkimage(i,:,3));
%             %   disp(['Error rate: ' num2str(error_rate)]);
%             % end
%             imwrite(uint8(extractwatermarkimage),'image/ext/extractwater_JPEG.bmp');
%             ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
%             berval=mean(ber);
%             result(count,ending)=ncval;
%             ending=ending+1;
%             result(count,ending)=berval;
%             %             time1(count,ending)=time1;
%             %             time2(count,ending)=time2;
%             ending=ending+1;
%         end
%         % % 椒盐噪声攻击测试，分层攻击图像
%         watermarkedimage = imread("image/WatermarkedImage/watermarked.bmp");
%         eding=1;
%         for Q=0.001:0.002:0.005
%             for level=1:3
%                 na=imnoise(watermarkedimage(:,:,level),'salt & pepper',Q);
%                 watermarklevel=watermarkimage(:,:,level);
%                 bitwater(level,:)=gainlevelwatermark(watermarklevel);
%                 Hd=double(na);   %temp.bmp 为压缩后的图片
%                 watermark=reshape(bitwater(level,:),128,64);
%                 %                 [nrw,ncw,vim]=size(watermark);
%                 %                 [nrf,ncf]=size(host(:,:,level));
%                 %                 rng(nrw);
%                 %                 pn_key_r=randperm(floor(nrf/4),nrw);
%                 %                 rand_r=pn_key_r;%%选取嵌入水印的位置
%                 %                 rng(ncw);
%                 %                 pn_key_c=randperm(floor(ncf/4),ncw);
%                 %                 rand_c=pn_key_c;
%                 %参数的设置
%                 %                 lambda=216; D=0.032;
%                 %                 delta=0.5*D;
%                 %                 psi=D;  %这个参数是啥？？
%                 %                 rho=1.5*D; eta=2*D;
%                 %                 blockSize=4;
%                 iter=1;
%                 lambda=216; D=0.032;
%                 delta=0.5*D;
%                 psi=D;
%                 rho=1.5*D; eta=2*D;
%                 blockSize=4;
%                 % NumOfColumn=256;
%                 % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
%                 for k=1:128
%                     for l=1:64
%                         DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
%                         DA_average=mean(DA(:));
%                         DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
%                         [ud,sd,vd]=svd(DA_lambda);
%                         rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
%                         if rc<0
%                             ud(:,1)=-ud(:,1);
%                             vd(:,1)=-vd(:,1);
%                         end
%                         g_hat=ud(2,1)-ud(3,1);
%                         if abs(g_hat)<=rho
%                             if g_hat>=0
%                                 Dw(k,l)='1';
%                             else
%                                 Dw(k,l)='0';
%                             end
%                         else
%                             if g_hat>rho
%                                 Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
%                             elseif g_hat<-rho
%                                 Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
%                             end
%                         end
%                     end
%                 end
%                 wa=watermark(:)';
%                 p=0;
%                 for i=1:128
%                     for j=1:64
%                         if Dw(i,j)~=watermark(i,j)
%                             p=p+1;
%                         end
%                     end
%                 end
%                 ber(level)=p/(8192);
%                 Dw2=reshape(Dw,1,8192);
%                 dd='';
%                 for i=1:8192
%                     dd=strcat(dd,num2str(Dw2(i)));
%                 end
%                 extractwatermarkimage(:,:,level)=gainlevelimage(dd);
%                 %                 watermark_d=reshape(Dw,[128,64]);
% 
%             end
%             % for i=1:32
%             %  error_rate = calculate_error_rate(extractwatermarkimage(i,:,3),watermarkimage(i,:,3));
%             %   disp(['Error rate: ' num2str(error_rate)]);
%             % end
%             imwrite(uint8(extractwatermarkimage),'image/ext/extractwater_JPEG.bmp');
%             ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
%             berval=mean(ber);
%             result(count,ending)=ncval;
%             ending=ending+1;
%             result(count,ending)=berval;
%             %             time1(count,ending)=time1;
%             %             time2(count,ending)=time2;
%             ending=ending+1;
%         end
%         % 缩放攻击测试，分层攻击图像
%         watermarkedimage = imread("image/WatermarkedImage/watermarked.bmp");
%         eding=1;
%         for Q=0.5:3.5:4
%             ra=imresize(watermarkedimage,Q);
%             ra=imresize(ra,[512 512]);
%             attacked=double(ra);
%             for level=1:3
%                 watermarklevel=watermarkimage(:,:,level);
%                 bitwater(level,:)=gainlevelwatermark(watermarklevel);
%                 Hd=attacked(:,:,level);   %temp.bmp 为压缩后的图片
%                 watermark=reshape(bitwater(level,:),128,64);
%                 %                 [nrw,ncw,vim]=size(watermark);
%                 %                 [nrf,ncf]=size(host(:,:,level));
%                 %                 rng(nrw);
%                 %                 pn_key_r=randperm(floor(nrf/4),nrw);
%                 %                 rand_r=pn_key_r;%%选取嵌入水印的位置
%                 %                 rng(ncw);
%                 %                 pn_key_c=randperm(floor(ncf/4),ncw);
%                 %                 rand_c=pn_key_c;
%                 %参数的设置
%                 %                 lambda=216; D=0.032;
%                 %                 delta=0.5*D;
%                 %                 psi=D;  %这个参数是啥？？
%                 %                 rho=1.5*D; eta=2*D;
%                 %                 blockSize=4;
%                 iter=1;
%                 lambda=216; D=0.032;
%                 delta=0.5*D;
%                 psi=D;
%                 rho=1.5*D; eta=2*D;
%                 blockSize=4;
%                 % NumOfColumn=256;
%                 % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
%                 for k=1:128
%                     for l=1:64
%                         DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
%                         DA_average=mean(DA(:));
%                         DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
%                         [ud,sd,vd]=svd(DA_lambda);
%                         rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
%                         if rc<0
%                             ud(:,1)=-ud(:,1);
%                             vd(:,1)=-vd(:,1);
%                         end
%                         g_hat=ud(2,1)-ud(3,1);
%                         if abs(g_hat)<=rho
%                             if g_hat>=0
%                                 Dw(k,l)='1';
%                             else
%                                 Dw(k,l)='0';
%                             end
%                         else
%                             if g_hat>rho
%                                 Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
%                             elseif g_hat<-rho
%                                 Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
%                             end
%                         end
%                     end
%                 end
%                 wa=watermark(:)';
%                 p=0;
%                 for i=1:128
%                     for j=1:64
%                         if Dw(i,j)~=watermark(i,j)
%                             p=p+1;
%                         end
%                     end
%                 end
%                 ber(level)=p/(8192);
%                 Dw2=reshape(Dw,1,8192);
%                 dd='';
%                 for i=1:8192
%                     dd=strcat(dd,num2str(Dw2(i)));
%                 end
%                 extractwatermarkimage(:,:,level)=gainlevelimage(dd);
%                 %                 watermark_d=reshape(Dw,[128,64]);
% 
%             end
%             % for i=1:32
%             %  error_rate = calculate_error_rate(extractwatermarkimage(i,:,3),watermarkimage(i,:,3));
%             %   disp(['Error rate: ' num2str(error_rate)]);
%             % end
%             imwrite(uint8(extractwatermarkimage),'image/ext/extractwater_JPEG.bmp');
%             ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
%             berval=mean(ber);
%             result(count,ending)=ncval;
%             ending=ending+1;
%             result(count,ending)=berval;
%             %             time1(count,ending)=time1;
%             %             time2(count,ending)=time2;
%             ending=ending+1;
%         end
%         % Translation攻击测试，分层攻击图像
%         watermarkedimage = imread("image/WatermarkedImage/watermarked.bmp");
%         eding=1;
%         for x=-10:10:20
%             attack=TranslatingAttack(watermarkedimage,x,x);%平移到[x,y]   floor((M/10)*5)=50%
%             attacked=TranslatingAttack(attack,-x,-x);
%             for level=1:3
%                 watermarklevel=watermarkimage(:,:,level);
%                 bitwater(level,:)=gainlevelwatermark(watermarklevel);
%                 Hd=double(attacked(:,:,level));   %temp.bmp 为压缩后的图片
%                 watermark=reshape(bitwater(level,:),128,64);
%                 %                 [nrw,ncw,vim]=size(watermark);
%                 %                 [nrf,ncf]=size(host(:,:,level));
%                 %                 rng(nrw);
%                 %                 pn_key_r=randperm(floor(nrf/4),nrw);
%                 %                 rand_r=pn_key_r;%%选取嵌入水印的位置
%                 %                 rng(ncw);
%                 %                 pn_key_c=randperm(floor(ncf/4),ncw);
%                 %                 rand_c=pn_key_c;
%                 %参数的设置
%                 %                 lambda=216; D=0.032;
%                 %                 delta=0.5*D;
%                 %                 psi=D;  %这个参数是啥？？
%                 %                 rho=1.5*D; eta=2*D;
%                 %                 blockSize=4;
%                 iter=1;
%                 lambda=216; D=0.032;
%                 delta=0.5*D;
%                 psi=D;
%                 rho=1.5*D; eta=2*D;
%                 blockSize=4;
%                 % NumOfColumn=256;
%                 % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
%                 for k=1:128
%                     for l=1:64
%                         DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
%                         DA_average=mean(DA(:));
%                         DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
%                         [ud,sd,vd]=svd(DA_lambda);
%                         rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
%                         if rc<0
%                             ud(:,1)=-ud(:,1);
%                             vd(:,1)=-vd(:,1);
%                         end
%                         g_hat=ud(2,1)-ud(3,1);
%                         if abs(g_hat)<=rho
%                             if g_hat>=0
%                                 Dw(k,l)='1';
%                             else
%                                 Dw(k,l)='0';
%                             end
%                         else
%                             if g_hat>rho
%                                 Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
%                             elseif g_hat<-rho
%                                 Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
%                             end
%                         end
%                     end
%                 end
%                 wa=watermark(:)';
%                 p=0;
%                 for i=1:128
%                     for j=1:64
%                         if Dw(i,j)~=watermark(i,j)
%                             p=p+1;
%                         end
%                     end
%                 end
%                 ber(level)=p/(8192);
%                 Dw2=reshape(Dw,1,8192);
%                 dd='';
%                 for i=1:8192
%                     dd=strcat(dd,num2str(Dw2(i)));
%                 end
%                 extractwatermarkimage(:,:,level)=gainlevelimage(dd);
%                 %                 watermark_d=reshape(Dw,[128,64]);
% 
%             end
%             % for i=1:32
%             %  error_rate = calculate_error_rate(extractwatermarkimage(i,:,3),watermarkimage(i,:,3));
%             %   disp(['Error rate: ' num2str(error_rate)]);
%             % end
%             imwrite(uint8(extractwatermarkimage),'image/ext/extractwater_JPEG.bmp');
%             ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
%             berval=mean(ber);
%             result(count,ending)=ncval;
%             ending=ending+1;
%             result(count,ending)=berval;
%             %             time1(count,ending)=time1;
%             %             time2(count,ending)=time2;
%             ending=ending+1;
%         end
%         % % Rotating攻击测试，分层攻击图像
%         watermarkedimage = imread("image/WatermarkedImage/watermarked.bmp");
%         eding=1;
%         for  Q=15:15:30
%             temp=RotatingAttack(watermarkedimage,Q);
% %             attacked=CorrectingGeometricAttacks(temp,512,512);
%             for level=1:3
%                 watermarklevel=watermarkimage(:,:,level);
%                 bitwater(level,:)=gainlevelwatermark(watermarklevel);
%                 Hd=double(temp(:,:,level));   %temp.bmp 为压缩后的图片
%                 watermark=reshape(bitwater(level,:),128,64);
%                 %                 [nrw,ncw,vim]=size(watermark);
%                 %                 [nrf,ncf]=size(host(:,:,level));
%                 %                 rng(nrw);
%                 %                 pn_key_r=randperm(floor(nrf/4),nrw);
%                 %                 rand_r=pn_key_r;%%选取嵌入水印的位置
%                 %                 rng(ncw);
%                 %                 pn_key_c=randperm(floor(ncf/4),ncw);
%                 %                 rand_c=pn_key_c;
%                 %参数的设置
%                 %                 lambda=216; D=0.032;
%                 %                 delta=0.5*D;
%                 %                 psi=D;  %这个参数是啥？？
%                 %                 rho=1.5*D; eta=2*D;
%                 %                 blockSize=4;
%                 iter=1;
%                 lambda=216; D=0.032;
%                 delta=0.5*D;
%                 psi=D;
%                 rho=1.5*D; eta=2*D;
%                 blockSize=4;
%                 % NumOfColumn=256;
%                 % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
%                 for k=1:128
%                     for l=1:64
%                         DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
%                         DA_average=mean(DA(:));
%                         DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
%                         [ud,sd,vd]=svd(DA_lambda);
%                         rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
%                         if rc<0
%                             ud(:,1)=-ud(:,1);
%                             vd(:,1)=-vd(:,1);
%                         end
%                         g_hat=ud(2,1)-ud(3,1);
%                         if abs(g_hat)<=rho
%                             if g_hat>=0
%                                 Dw(k,l)='1';
%                             else
%                                 Dw(k,l)='0';
%                             end
%                         else
%                             if g_hat>rho
%                                 Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
%                             elseif g_hat<-rho
%                                 Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
%                             end
%                         end
%                     end
%                 end
%                 wa=watermark(:)';
%                 p=0;
%                 for i=1:128
%                     for j=1:64
%                         if Dw(i,j)~=watermark(i,j)
%                             p=p+1;
%                         end
%                     end
%                 end
%                 ber(level)=p/(8192);
%                 Dw2=reshape(Dw,1,8192);
%                 dd='';
%                 for i=1:8192
%                     dd=strcat(dd,num2str(Dw2(i)));
%                 end
%                 extractwatermarkimage(:,:,level)=gainlevelimage(dd);
%                 %                 watermark_d=reshape(Dw,[128,64]);
% 
%             end
%             % for i=1:32
%             %  error_rate = calculate_error_rate(extractwatermarkimage(i,:,3),watermarkimage(i,:,3));
%             %   disp(['Error rate: ' num2str(error_rate)]);
%             % end
%             imwrite(uint8(extractwatermarkimage),'image/ext/extractwater_JPEG.bmp');
%             ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
%             berval=mean(ber);
%             result(count,ending)=ncval;
%             ending=ending+1;
%             result(count,ending)=berval;
%             %             time1(count,ending)=time1;
%             %             time2(count,ending)=time2;
%             ending=ending+1;
%         end
%         % % 剪切攻击测试
%         watermarkedimage = imread("image/WatermarkedImage/watermarked.bmp");
%         eding=1;
%         for  Q=1:1:2
%             attacked=watermarkedimage;
%             attacked(1:64*Q,:,:)=255;
%             for level=1:3
%                 watermarklevel=watermarkimage(:,:,level);
%                 bitwater(level,:)=gainlevelwatermark(watermarklevel);
%                 Hd=double(attacked(:,:,level));   %temp.bmp 为压缩后的图片
%                 watermark=reshape(bitwater(level,:),128,64);
%                 %                 [nrw,ncw,vim]=size(watermark);
%                 %                 [nrf,ncf]=size(host(:,:,level));
%                 %                 rng(nrw);
%                 %                 pn_key_r=randperm(floor(nrf/4),nrw);
%                 %                 rand_r=pn_key_r;%%选取嵌入水印的位置
%                 %                 rng(ncw);
%                 %                 pn_key_c=randperm(floor(ncf/4),ncw);
%                 %                 rand_c=pn_key_c;
%                 %参数的设置
%                 %                 lambda=216; D=0.032;
%                 %                 delta=0.5*D;
%                 %                 psi=D;  %这个参数是啥？？
%                 %                 rho=1.5*D; eta=2*D;
%                 %                 blockSize=4;
%                 iter=1;
%                 lambda=216; D=0.032;
%                 delta=0.5*D;
%                 psi=D;
%                 rho=1.5*D; eta=2*D;
%                 blockSize=4;
%                 % NumOfColumn=256;
%                 % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
%                 for k=1:128
%                     for l=1:64
%                         DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
%                         DA_average=mean(DA(:));
%                         DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
%                         [ud,sd,vd]=svd(DA_lambda);
%                         rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
%                         if rc<0
%                             ud(:,1)=-ud(:,1);
%                             vd(:,1)=-vd(:,1);
%                         end
%                         g_hat=ud(2,1)-ud(3,1);
%                         if abs(g_hat)<=rho
%                             if g_hat>=0
%                                 Dw(k,l)='1';
%                             else
%                                 Dw(k,l)='0';
%                             end
%                         else
%                             if g_hat>rho
%                                 Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
%                             elseif g_hat<-rho
%                                 Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
%                             end
%                         end
%                     end
%                 end
%                 wa=watermark(:)';
%                 p=0;
%                 for i=1:128
%                     for j=1:64
%                         if Dw(i,j)~=watermark(i,j)
%                             p=p+1;
%                         end
%                     end
%                 end
%                 ber(level)=p/(8192);
%                 Dw2=reshape(Dw,1,8192);
%                 dd='';
%                 for i=1:8192
%                     dd=strcat(dd,num2str(Dw2(i)));
%                 end
%                 extractwatermarkimage(:,:,level)=gainlevelimage(dd);
%                 %                 watermark_d=reshape(Dw,[128,64]);
% 
%             end
%             % for i=1:32
%             %  error_rate = calculate_error_rate(extractwatermarkimage(i,:,3),watermarkimage(i,:,3));
%             %   disp(['Error rate: ' num2str(error_rate)]);
%             % end
%             imwrite(uint8(extractwatermarkimage),'image/ext/extractwater_JPEG.bmp');
%             ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
%             berval=mean(ber);
%             result(count,ending)=ncval;
%             ending=ending+1;
%             result(count,ending)=berval;
%             %             time1(count,ending)=time1;
%             %             time2(count,ending)=time2;
%             ending=ending+1;
%         end
%         count=count+1;
%     end
%     %     result(count,1)=TT;
%     %     result(count,2)=mean(PSNR);
%     %     result(count,3)=mean(SSIM);
%     %     result(count,4)=mean(NC);
%     %     result(count,5)=mean(berval);
%     % result(count,6)=mean(time1);
%     % result(count,7)=mean(time2);
%     %     A=[mean(AttackNC1),mean(AttackNC2),mean(AttackNC3),mean(AttackNC4),mean(AttackNC5),mean(AttackNC6),mean(AttackNC7),mean(AttackNC8),mean(AttackNC9)];
%     %     result(count,6)=mean(A);
%     %     count=count+1;
%     %     clear AttackNC1;
%     %     clear AttackNC2;
%     %     clear AttackNC3;
%     %     clear AttackNC4;
%     %     clear AttackNC5;
%     %     clear AttackNC6;
%     %     clear AttackNC7;
%     %     clear AttackNC8;
%     %     clear AttackNC9;
    end
end
xlswrite("mydataCompare.xlsx", result, "sheet5");

        watermarkedimage = imread("image/WatermarkedImage/watermarked.bmp");
        eding=1;
        for  Q=1:1:2
            attacked=watermarkedimage;
            attacked(1:64*Q,:,:)=0;
            for level=1:3
                watermarklevel=watermarkimage(:,:,level);
                bitwater(level,:)=gainlevelwatermark(watermarklevel);
                Hd=double(attacked(:,:,level));   %temp.bmp 为压缩后的图片
                watermark=reshape(bitwater(level,:),128,64);
                %                 [nrw,ncw,vim]=size(watermark);
                %                 [nrf,ncf]=size(host(:,:,level));
                %                 rng(nrw);
                %                 pn_key_r=randperm(floor(nrf/4),nrw);
                %                 rand_r=pn_key_r;%%选取嵌入水印的位置
                %                 rng(ncw);
                %                 pn_key_c=randperm(floor(ncf/4),ncw);
                %                 rand_c=pn_key_c;
                %参数的设置
                %                 lambda=216; D=0.032;
                %                 delta=0.5*D;
                %                 psi=D;  %这个参数是啥？？
                %                 rho=1.5*D; eta=2*D;
                %                 blockSize=4;
                iter=1;
                lambda=216; D=0.032;
                delta=0.5*D;
                psi=D;
                rho=1.5*D; eta=2*D;
                blockSize=4;
                % NumOfColumn=256;
                % Hd=u_o(:,1:NumOfColumn)*s_o(1:NumOfColumn,1:NumOfColumn)*v_o(:,1:NumOfColumn)';
                for k=1:128
                    for l=1:64
                        DA=Hd((k-1)*blockSize+1:(k-1)*blockSize+blockSize,(l-1)*blockSize+1:(l-1)*blockSize+blockSize);
                        DA_average=mean(DA(:));
                        DA_lambda=DA+(lambda-DA_average)*ones(blockSize)+(rand(blockSize)-0.5);
                        [ud,sd,vd]=svd(DA_lambda);
                        rc=ud(:,1)'*DA_lambda*ones(blockSize,1)+ones(1,blockSize)*DA_lambda*vd(:,1);
                        if rc<0
                            ud(:,1)=-ud(:,1);
                            vd(:,1)=-vd(:,1);
                        end
                        g_hat=ud(2,1)-ud(3,1);
                        if abs(g_hat)<=rho
                            if g_hat>=0
                                Dw(k,l)='1';
                            else
                                Dw(k,l)='0';
                            end
                        else
                            if g_hat>rho
                                Dw(k,l)=num2str(mod(floor((g_hat-eta)/D+0.5), 2));
                            elseif g_hat<-rho
                                Dw(k,l)=num2str(mod(ceil((g_hat+eta)/D-0.5)-1, 2));
                            end
                        end
                    end
                end
                wa=watermark(:)';
                p=0;
                for i=1:128
                    for j=1:64
                        if Dw(i,j)~=watermark(i,j)
                            p=p+1;
                        end
                    end
                end
                ber(level)=p/(8192);
                Dw2=reshape(Dw,1,8192);
                dd='';
                for i=1:8192
                    dd=strcat(dd,num2str(Dw2(i)));
                end
                extractwatermarkimage(:,:,level)=gainlevelimage(dd);
                %                 watermark_d=reshape(Dw,[128,64]);
                
            end
            % for i=1:32
            %  error_rate = calculate_error_rate(extractwatermarkimage(i,:,3),watermarkimage(i,:,3));
            %   disp(['Error rate: ' num2str(error_rate)]);
            % end
            imwrite(uint8(extractwatermarkimage),'image/ext/extractwater_crop.bmp');
            ncval=colornc(uint8(extractwatermarkimage),uint8(watermarkimage));
            berval=mean(ber);
            result(count,ending)=ncval;
            ending=ending+1;
            result(count,ending)=berval;
            figure();
            subplot(122),imshow('image/ext/extractwater_crop.bmp'),title([' crop NC/BER=',num2str(ncval),'/',num2str(berval),]);
            %             time1(count,ending)=time1;
            %             time2(count,ending)=time2;
            ending=ending+1;
        end