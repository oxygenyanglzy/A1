clc              %Proposed
clear all;
blocksize=4;
load position;
a=1;b=-1;c=-1;d=0;
key=6;
im={'barnfall.jpg'};%{'f16.jpg','safari10.jpg','woman.jpg','pine.jpg','baboon.jpg','butrfly1.jpg','house.jpg','barnfall.jpg'};%,};
% for imm=1:length(im)
k=1;%
for imm=1:length(im)
    % k=1;
    %m=imm;
    m=1;
    for Qfac=26%2:2:34   %Qfac为量化步长
        % T=0.012
        %实验结果表明：block=4×4具有很强的不可见性，
        %优越性明显；block=8×8具有较好的鲁棒性，效果与前者相差无几，
        %在水印容量较大时可以考虑用较小的分块尺寸得到多个分块用来隐藏水印。
        %改变分别改变U和
        %H='lena.jpg';%;'peppers.jpg';%'baboon.jpg';%'F16.jpg';%
        %  t1=0;
        %  t2=0;
        % t1=cputime;
        H=im{1,imm};
        %H='lena.jpg';
        %W='3073232.tif';
        W1='20010073260.0.jpg';
        W='200100732.jpg';
        %W='colorwatermark32.bmp';
        %W='beijing1.bmp';
        %W='ldu3232.jpg';
        Host=imread(H);
        Water=imread(W1);
        Water1=imread(W);
        lenw=size(Water,1)*size(Water,2)*8;
        figure(1),subplot(121),imshow(Host),title('Host image');
        subplot(122),imshow(Water),title('watermark image');
        %水印Affine变换
        %--------------------------------------------------
        %原始水印处理 调用子函数形成分层水印levelwatermark
        %water307=imread('colorwatermark32.bmp');
        %tic;
        for level=1:3
            %  levelwatermark(level,:)=gainlevelwatermark(Water(:,:,level));
            temp1=Water(:,:,level);
            afterAffine=Affine(temp1,a,b,c,d,key,0);
            levelwatermark(level,:)=gainlevelwatermark(afterAffine);
        end
        %--------------------------------------------------
        %水印嵌入
        %psnrval=0;
        
        
        for level=1:3
            % level=3
            Hostlevel=double(Host(:,:,level));
            Hostlevelback=Hostlevel;
            nlevels =1;%[0, 1, 3] ;        % Decomposition level
            pfilter = 'maxflat' ;%'9-7';%              % Pyramidal filter
            dfilter = 'dmaxflat7' ; %  'pkva';           % Directional filter
            
            % Nonsubsampled Contourlet decomposition
            coeffs = nsctdec( Hostlevel, nlevels, dfilter, pfilter );
            
            Hostlevel=coeffs{1,1};%%save in su1
            
            Waterold=levelwatermark(level,:);
            [Hr,Hc]=size(Hostlevel);
            waterposition=1;
            
            for  counter=1:lenw
                i=position(1,counter);%根据position取出应该的嵌入位置
                j=position(2,counter);
                watermark=Waterold(1,counter);%%根据水印序列取出水印位
                block=Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
                %将经过变换后的
                [u,s]=schur(block);
                lmt= round(s(1,1)/Qfac);
                if watermark=='1'
                    wat=1;
                else
                    wat=0;
                end
                Econd=xor(mod(lmt,2),wat);
                if Econd==0
                    s(1,1)=lmt*Qfac+(Qfac/2);
                else
                    s(1,1)=lmt*Qfac-(Qfac/2);
                end
                blocknew=u*s*u';
                Hostlevel((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize)=blocknew;
            end
            % Reconstruct image
            coeffs{1,1}=Hostlevel;%%su2
            Hostlevel= nsctrec(coeffs, dfilter, pfilter) ;
            watermarkedim(:,:,level)=Hostlevel;
            imwrite(uint8(Hostlevel),'www.bmp');
            imwrite(Host(:,:,level),'hhh.bmp');
        end
        
        imwrite(uint8(watermarkedim),'watermarked.bmp');
        psnrval(k,m)=colorpsnr(watermarkedim,Host);
        ssimval(k,m)=colorssim(watermarkedim,Host);
        hand=figure(2);
        subplot(431),imshow(Host),title('Original image');
        subplot(432),imshow('watermarked.bmp'),title([' PSNR=',num2str(psnrval(k,m)),' SSIM=',num2str(ssimval(k,m)),' Qfac=',num2str(Qfac)]);
        %toc;
        % t2=cputime;
        
        %%%%--------------------------------------------------
        %%水印提取
        %tic;
        watermarkedim=double(imread('watermarked.bmp'));
        for level=1:3
            wa=watermarkedim(:,:,level);
            Waterold=levelwatermark(level,:);
            % Nonsubsampled Contourlet decomposition
            coeffs = nsctdec( wa, nlevels, dfilter, pfilter );
            wa=coeffs{1,1};%%su3
            for counter=1:lenw
                i=position(1,counter);
                j=position(2,counter);
                block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
                [u1,s1]=schur(block);
                Rext=fix(s1(1,1)/Qfac);
                ext=mod(Rext,2);
                if ext==0
                    watermark='0';
                else
                    watermark='1';
                end
                ExWater(1,counter)=watermark;
            end
            wateronedim(level,:)=ExWater;
        end
        %--------------------------------------------------
        %水印恢复
        %extractlevelwatermark=watermarkrestore(wateronedim);
        for level=1:3
            extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
            extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
        end
        %toc;
        imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
        ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water1));
        subplot(133),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
        % name=[H(1:length(H)-4),'Unattack',W(1:length(W)-4),'.fig'];
        % saveas(hand,name);
        % extractt(k,m)=cputime-t2;
        
        % %% JPEG2000 测试 一  图像整体攻击提取效果较差
        %  watermarkedim=imread('watermarked.bmp');
        %  for QA=4;%1:3:7%:10
        %   for level=1:3
        %       imwrite(watermarkedim(:,:,level),'temp.bmp','jp2','compressionratio',QA);
        %       wb=double(imread('temp.bmp'));
        %       ttttt(:,:,level)=wb;
        % %Nonsubsampled Contourlet decomposition
        %     coeffs = nsctdec( wb, nlevels, dfilter, pfilter );
        %      wa=coeffs{1,1};%%su3
        %     for  counter=1:lenw
        %         i=position(1,counter);
        %         j=position(2,counter);
        %         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        %          [u2,s2]=schur(block);
        %         Rext=fix(s2(1,1)/Qfac);
        %         ext=mod(Rext,2);
        %         if ext==0
        %             watermark='0';
        %         else
        %             watermark='1';
        %         end
        %         ExWater(1,counter)=watermark;
        %
        %      end
        %      wateronedim(level,:)=ExWater;
        %
        % end
        % %
        % % %水印回复
        % for level=1:3
        %         extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
        %         extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
        %  end
        % imwrite(ttttt,'temp.bmp');
        % psnrval(k,m)=colorpsnr(ttttt,Host);
        % imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
        % ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
        % hand=figure(QA+2);
        % subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
        % subplot(122),imshow(uint8(ttttt)),title([num2str(QA),'  PSNR=',num2str(psnrval(k,m))]);
        % name=[H(1:length(H)-4),'JPEG2000_', num2str(QA),'-',W(1:length(W)-4),'.fig'];
        % % saveas(hand,name);
        % % result(1,QA)=psnrval;
        % % result(2,QA)=ncval;
        % m=m+1;
        % end
        % %
        
        
        
        % %
        % % %
        % % %---------------------------------------------------------------------
        %JPEG 测试 二  图像分层攻击提取效果较好
        
         watermarkedim=imread('watermarked.bmp');
         for QA=70%10:30:70
            for level=1:3
              imwrite(watermarkedim(:,:,level),'temp.bmp','jpg','quality',QA);
              wb=double(imread('temp.bmp'));
              ttttt(:,:,level)=wb;
        
        % Nonsubsampled Contourlet decomposition
             coeffs = nsctdec( wb, nlevels, dfilter, pfilter );
             wa=coeffs{1,1};%%su3
            for  counter=1:lenw
                 i=position(1,counter);
                 j=position(2,counter);
                 block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
                   [u2,s2]=schur(block);
                Rext=fix(s2(1,1)/Qfac);
                ext=mod(Rext,2);
                if ext==0
                    watermark='0';
                else
                    watermark='1';
                end
                ExWater(1,counter)=watermark;
            end
           wateronedim(level,:)=ExWater;
        
           end
        %
        %水印回复
        
        for level=1:3
                extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
                extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
         end
        
        imwrite(ttttt,'temp.bmp');
        psnrval(k,m)=colorpsnr(ttttt,Host);
        imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
        ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water1));
        hand=figure(QA/10+2);
        subplot(121),imshow( uint8(extractlevelwatermark)),title(['Delamination Extracted Watermark nc=',num2str(ncval(k,m))]);
        subplot(122),imshow(uint8(ttttt)),title([num2str(QA),'  PSNRr=',num2str(psnrval(k,m))]);
        name=[H(1:length(H)-4),'JPEG_', num2str(QA),'-',W(1:length(W)-4),'.fig'];
        m=m+1;
        % saveas(hand,name);
        % result(1,QA/10)=psnrval;
        % result(2,QA/10)=ncval;
        end
        
        
        %椒盐噪声&&高斯噪声
        a=imread('watermarked.bmp');
        for QQ=2%1:1:3
            %wb=imnoise(a,'salt & pepper',0.002*QQ);  %salt & pepper
            wb=imnoise(a,'gaussian',0,0.001*QQ); %gausian
           imwrite(wb,'temp.bmp');
           imshow(wb);
          for level=1:3
             ws=double(wb(:,:,level));
             coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
             wa=coeffs{1,1};%%su3
        for  counter=1:lenw
            i=position(1,counter);
            j=position(2,counter);
        
                block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
                [u2,s2]=schur(block);
                Rext=fix(s2(1,1)/Qfac);
                ext=mod(Rext,2);
                if ext==0
                    watermark='0';
                else
                    watermark='1';
                end
                ExWater(1,counter)=watermark;
             end
             wateronedim(level,:)=ExWater;
        
        end
        % %
        %水印回复
        
        for level=1:3
                extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
                extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
         end
        
        psnrval(k,m)=colorpsnr(wb,Host);
        imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
        ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water1));
        hand=figure(333);
        subplot(121),imshow( uint8(extractlevelwatermark)),title(['NOISE Extracted Watermark nc=',num2str(ncval(k,m))]);
        subplot(122),imshow(uint8(wb)),title([num2str(QQ),' PSNR=',num2str(psnrval(k,m))]);
        name=[H(1:length(H)-4),'salt_pepper-', num2str(0.01*QQ),'-',W(1:length(W)-4),'.fig'];
        % saveas(hand,name);
        % result(1,QQ)=psnrval;
        % result(2,QQ)=ncval;
        m=m+1;
        end
        
%         %%%%%%%%%%%%%%%%%%
%         watermarkedim=imread('watermarked.bmp');
%          for QQ=3%1:2:5
%         
%           %h = fspecial('average');%('gaussian');%, 50, 45);%%('average')%,'disk','gaussian',
%         
%           for level=1:3
%               %%wb = imfilter(watermarkedim(:,:,level), h);
%               wb=medfilt2(watermarkedim(:,:,level),[QQ,1]);
%                ttttt(:,:,level)=wb;
%             ws=double(wb);
%             waterposition=1;
%              coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
%              wa=coeffs{1,1};%%su3
%           for  counter=1:lenw
%             i=position(1,counter);
%             j=position(2,counter);
%         
%                 block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
%                 [u2,s2]=schur(block);
%                 Rext=fix(s2(1,1)/Qfac);
%                 ext=mod(Rext,2);
%                 if ext==0
%                     watermark='0';
%                 else
%                     watermark='1';
%                 end
%                 ExWater(1,counter)=watermark;
%            end
%              wateronedim(level,:)=ExWater;
%         
%            end
%         
%         %%%水印回复
%         
%          for level=1:3
%                 extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
%                 extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
%          end
%         imwrite(ttttt,'temp.bmp');
%         psnrval(k,m)=colorpsnr(ttttt,Host);
%         imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
%         ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
%         hand=figure(QQ+2);
%         subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
%         subplot(122),imshow(uint8(ttttt)),title([num2str(QQ),'  PSNRr=',num2str(psnrval(k,m))]);
%         name=[H(1:length(H)-4),'Median-', num2str(QQ),'-',W(1:length(W)-4),'.fig'];
%         % saveas(hand,name);
%         % result(1,QQ)=psnrval;
%         % result(2,QQ)=ncval;
%         m=m+1;
%          end
        %%% % %%%%% Butterlow-pass filtering
        %
        % watermarkedim=imread('watermarked.bmp');
        %  for QQ=5%1:5%1:5
        %     for level=1:3
        %         Ib=watermarkedim(:,:,level);
        %         f=double(Ib);
        %         g=fft2(f);
        %         g=fftshift(g);
        %         [N1,N2]=size(g);
        %         d0=100;
        %         n1=fix(N1/2);
        %         n2=fix(N2/2);
        %         for i=1:N1
        %             for j=1:N2
        %                 d=sqrt((i-n1)^2+(j-n2)^2);
        %                 h=1/(1+0.414*(d/d0)^(2*QQ));
        %                 result(i,j)=h*g(i,j);
        %             end
        %         end
        %         result=ifftshift(result);
        %         X2=ifft2(result);
        %         X3=uint8(real(X2));
        %         ttttt(:,:,level)=X3;
        %         wb=double(X3);
        %             % Nonsubsampled Contourlet decomposition
        %     coeffs = nsctdec( wb, nlevels, dfilter, pfilter );
        %     wa=coeffs{1,1};%%su3
        %     for  counter=1:lenw
        %     i=position(1,counter);
        %     j=position(2,counter);
        %         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        %         [u2,s2]=schur(block);
        %         Rext=fix(s2(1,1)/Qfac);
        %         ext=mod(Rext,2);
        %         if ext==0
        %             watermark='0';
        %         else
        %             watermark='1';
        %         end
        %         ExWater(1,counter)=watermark;
        %     end
        %          wateronedim(level,:)=ExWater;
        %     end
        % %                  p=0;
        % %                 for level=1:3
        % %                     for i=1:8192
        % %                         if wateronedim(level,i)~=levelwatermark(level,i)
        % %                             p=p+1;
        % %                         end
        % %                     end
        % %                 end
        % %                 ber(k,m)=p/24576;
        % % %                % disp(ber(m,k));
        % % %                 k=k+1;
        % % %%%水印回复
        % %
        % % %%水印回复
        %   for level=1:3
        %         extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
        %         extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
        %  end
        % imwrite(ttttt,'temp.bmp');
        % psnrval(k,m)=colorpsnr(ttttt,Host);
        % imwrite(uint8(extractlevelwatermark),'extrwater.bmp')
        % % name1=[H(1:length(H)-4),'ButterLowpass-', num2str(QQ),'-',W(1:length(W)-4),'.bmp'];
        % % imwrite(uint8(extractlevelwatermark),['D:/attackimage/md5schur',name1]);
        % ncval(k,m)=colornc(uint8(extractlevelwatermark),uint8(Water));
        % hand=figure(QQ);
        % subplot(121),imshow( uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
        % subplot(122),imshow(uint8(ttttt)),title([num2str(QQ),'  PSNRr=',num2str(psnrval(k,m))]);
        % name=[H(1:length(H)-4),'ButterLowpass', '-',num2str(QQ),'-',W(1:length(W)-4),'.fig']; saveas(hand,name);
        %  m=m+1;
        %  end
        % a=imread('watermarked.bmp');
        % for QQ=8%2:2:12
        %    wb=imresize(a,QQ*0.5);
        %    wb=imresize(wb,[512 512]);
        %    imwrite(wb,'temp.bmp');
        %    for level=1:3
        %     ws=double(wb(:,:,level));
        %     waterposition=1;
        %      coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
        %      wa=coeffs{1,1};
        %  for  counter=1:lenw
        %     i=position(1,counter);
        %     j=position(2,counter);
        %         block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
        %         [u2,s2]=schur(block);
        %         Rext=fix(s2(1,1)/Qfac);
        %         ext=mod(Rext,2);
        %         if ext==0
        %             watermark='0';
        %         else
        %             watermark='1';
        %         end
        %         ExWater(1,counter)=watermark;
        %   end
        %        wateronedim(level,:)=ExWater;
        %
        %    end
        %
        % %
        % %        p=0;
        % %                 for level=1:3
        % %                     for i=1:8192
        % %                         if wateronedim(level,i)~=levelwatermark(level,i)
        % %                             p=p+1;
        % %                         end
        % %                     end
        % %                 end
        % %                 ber(k,m)=p/24576;
        % % %                % disp(ber(m,k));
        % % %                % k=k+1;
        % %水印回复
        % for level=1:3
        %         extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
        %         extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
        % end
        % psnrval(k,m)=colorpsnr(wb,Host);
        % imwrite(uint8(extractlevelwatermark),'extrwater.jpg')
        % ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
        % hand=figure(QQ+2);
        % subplot(121),imshow(uint8(extractlevelwatermark)),title(['Extracted Watermark nc=',num2str(ncval(k,m))]);
        % subplot(122),imshow(uint8(wb)),title([num2str(0.5*QQ),'  PSNRr=',num2str(psnrval(k,m))]);
        % name=[H(1:length(H)-4),'Resize_', num2str(QQ),'-',W(1:length(W)-4),'.fig'];
        % m=m+1;
        %  end
        
        
        %%剪切攻击
        for QQ=2%1:4
            wb=imread('watermarked.bmp');
            %    wb(1:64*QQ,1:64*QQ,:)=0; % Cropping attack
            wb(:,1:64*QQ,:)=0;
            imwrite(wb,'temp.bmp');
            
            for level=1:3
                ws=double(wb(:,:,level));
                coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
                wa=coeffs{1,1};%%su3
                for  counter=1:lenw
                    i=position(1,counter);
                    j=position(2,counter);
                    
                    block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
                    [u2,s2]=schur(block);
                    Rext=fix(s2(1,1)/Qfac);
                    ext=mod(Rext,2);
                    if ext==0
                        watermark='0';
                    else
                        watermark='1';
                    end
                    ExWater(1,counter)=watermark;
                end
                wateronedim(level,:)=ExWater;
                
            end
            %        p=0;
            %                 for level=1:3
            %                     for i=1:8192
            %                         if wateronedim(level,i)~=levelwatermark(level,i)
            %                             p=p+1;
            %                         end
            %                     end
            %                 end
            %                 ber(k,m)=p/24576;
            %                % disp(ber(m,k));
            
            %水印回复
            for level=1:3
                extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
                extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
            end
            %extractlevelwatermark=watermarkrestore(wateronedim);
            psnrval(k,m)=colorpsnr(wb,Host);
            imwrite(uint8(extractlevelwatermark),'extrwater.bmp');
            % name1=[H(1:length(H)-4),'Cropping__', num2str(QQ),'-',W(1:length(W)-4),'.bmp'];
            %  imwrite(uint8(extractlevelwatermark),['D:/attackimage/md5schur',name1]);
            ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
            % hand=figure(1+QQ);
            figure(33),subplot(121),imshow( uint8(extractlevelwatermark)),title(['Cropping Extracted Watermark nc=',num2str(ncval(k,m))]);
            figure(33),subplot(122),imshow(uint8(wb)),title([num2str(QQ),'  PSNRr=',num2str(psnrval(k,m))]);
            name=[H(1:length(H)-4),'Crop_', num2str(QQ*0.125*0.5),'-',W(1:length(W)-4),'.fig'];
            saveas(hand,name);
        end
        
        
        load geoData512
        load position
        watermarkedimg=imread('watermarked.bmp');
        %---------------几何攻击_旋转和矫正---------------------
        % % 512×512+消锯齿
        % watermarkedimg=imread('f1646-num.bmp');
        tic
        %[attackimg,aname]=GeometricAttack(watermarkedimg);
        attackimg=RotatingAttack(watermarkedimg,30);
        figure(3),subplot(341),imshow(attackimg),title('RotatingAttack image');
        % % 找到四个点
        watermarked=Imagecorrection(attackimg);
        for level=1:3
            ws=double(watermarked(:,:,level));
            coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
            wa=coeffs{1,1};%%su3
            for  counter=1:lenw
                i=position(1,counter);
                j=position(2,counter);
                block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
                [u2,s2]=schur(block);
                Rext=fix(s2(1,1)/Qfac);
                ext=mod(Rext,2);
                if ext==0
                    watermark='0';
                else
                    watermark='1';
                end
                ExWater(1,counter)=watermark;
            end
            wateronedim(level,:)=ExWater;
            
        end
        %水印恢复
        %extractlevelwatermark=watermarkrestore(wateronedim);
        for level=1:3
            extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
            extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
        end
        %toc;
        imwrite(uint8(extractlevelwatermark),'RotatingAttack extrwater.bmp')
        ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
        figure(3),subplot(342),imshow(uint8(watermarked)),title('RotatingAttackcorrected image2');
        figure(4),subplot(221),imshow( uint8(extractlevelwatermark)),title(['RotatingAttack Extracted Watermark nc=',num2str(ncval(k,m))]);
%        [cornerpoint,attacktype]=LocatingVertices4(attackimg);
%         if attacktype %平移
%             watermarked=CorrectingTranslate2(attackimg,cornerpoint);
%         else %其他情况
%             src1=[cornerpoint(:,:,1);cornerpoint(:,:,2);cornerpoint(:,:,3);cornerpoint(:,:,4)];%改，根据不同角度确定对应的是cornerpoint的哪个点
%             pt_M=512;
%             tar1=[1,1;pt_M,1;pt_M,pt_M;1,pt_M];
%             TForm1 = fitgeotrans(src1,tar1,'Projective'); %Projective
%             attackimg_pt = imwarp(attackimg, TForm1);
%             %     imwrite(attackimg_pt,'attackimg_pt.bmp');
%             %     imshow(attackimg_pt);
%             % % 去黑边及消除锯齿效应
%             watermarked=RemovingBlackEdges4(attackimg_pt,512,512); % 去黑边
%             %     imwrite(uint8(watermarked),'removeblackedge_F16.bmp');
%             figure(3),subplot(122),imshow(watermarked),title('remove black matte image');
%             watermarked=sawtoothProcess(double(watermarked)); %消除锯齿效应
%             %     imwrite(uint8(watermarked),'reducesawtootheffect.bmp');
%             figure(4),subplot(121),imshow(uint8(watermarked)),title('remove sawtooth image');
%         end
        geotime=toc;
        % attname=['./geo_graph/',aname,'.bmp']; %/geo_graph/rotate30.bmp'
        % corrname=['./geo_graph/cor-',aname,'.bmp']; %/geo_graph/cor-rotate30.bmp'
        % imwrite(uint8(attackimg),attname);
        % imwrite(uint8(watermarked),corrname);
        % psnrval=colorpsnr(uint8(watermarked),watermarkedimg);
        % geoData512{end+1,1}=aname;
        % geoData512{end,2}=psnrval;
        % geoData512{end,5}=geotime;
        
        %------------------------------------------------------
        
        %---------------几何攻击_倾斜和矫正---------------------
        tar1=[50,1;512,60;450,512;1,400];
        attackimg=SlantAttack(watermarkedimg,tar1);
        figure(3),subplot(343),imshow(attackimg),title('SlantAttack image');
        % % 找到四个点
        watermarked=Imagecorrection(attackimg);
        for level=1:3
            ws=double(watermarked(:,:,level));
            coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
            wa=coeffs{1,1};%%su3
            for  counter=1:lenw
                i=position(1,counter);
                j=position(2,counter);
                block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
                [u2,s2]=schur(block);
                Rext=fix(s2(1,1)/Qfac);
                ext=mod(Rext,2);
                if ext==0
                    watermark='0';
                else
                    watermark='1';
                end
                ExWater(1,counter)=watermark;
            end
            wateronedim(level,:)=ExWater;
            
        end
        %水印恢复
        %extractlevelwatermark=watermarkrestore(wateronedim);
        for level=1:3
            extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
            extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
        end
        %toc;
        imwrite(uint8(extractlevelwatermark),'SlantAttack extrwater.bmp')
        ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
        figure(3),subplot(344),imshow(uint8(watermarked)),title('SlantAttackcorrected image2');
        figure(4),subplot(222),imshow( uint8(extractlevelwatermark)),title(['SlantAttack Extracted Watermark nc=',num2str(ncval(k,m))]);
        %------------------------------------------------------
        
        %---------------几何攻击_仿射和矫正---------------------
        load position
        para=0.07;
        attackimg=AffineAttack(watermarkedimg,1,para,para);
        imwrite(uint8(attackimg),'AffineAttack.bmp')
        % % 找到四个点
        watermarked=Imagecorrection(attackimg);
        for level=1:3
            ws=double(watermarked(:,:,level));
            coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
            wa=coeffs{1,1};%%su3
            for  counter=1:lenw
                i=position(1,counter);
                j=position(2,counter);
                block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
                [u2,s2]=schur(block);
                Rext=fix(s2(1,1)/Qfac);
                ext=mod(Rext,2);
                if ext==0
                    watermark='0';
                else
                    watermark='1';
                end
                ExWater(1,counter)=watermark;
            end
            wateronedim(level,:)=ExWater;
            
        end
        %水印恢复
        %extractlevelwatermark=watermarkrestore(wateronedim);
        for level=1:3
            extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
            extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
        end
        %toc;
        imwrite(uint8(extractlevelwatermark),'AffineAttack extrwater.bmp')
        ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
        figure(3),subplot(345),imshow(attackimg),title('AffineAttack image');
        figure(3),subplot(346),imshow(uint8(watermarked)),title('AffineAttackcorrected image');
        imwrite(uint8(watermarked),'AffineAttackcorrected.bmp')
        figure(4),subplot(223),imshow( uint8(extractlevelwatermark)),title(['AffineAttack Extracted Watermark nc=',num2str(ncval(k,m))]);
        %------------------------------------------------------
        
        %---------------几何攻击_平移和矫正---------------------
        attackimg=TranslatingAttack(watermarkedimg,35,25);
        figure(3),subplot(347),imshow(attackimg),title('TranslatingAttack image');
        % % 找到四个点
        watermarked=Imagecorrection(attackimg);
        for level=1:3
            ws=double(watermarked(:,:,level));
            coeffs = nsctdec( ws, nlevels, dfilter, pfilter );
            wa=coeffs{1,1};%%su3
            for  counter=1:lenw
                i=position(1,counter);
                j=position(2,counter);
                block=wa((i-1)*blocksize+1:(i-1)*blocksize+blocksize,(j-1)*blocksize+1:(j-1)*blocksize+blocksize);
                [u2,s2]=schur(block);
                Rext=fix(s2(1,1)/Qfac);
                ext=mod(Rext,2);
                if ext==0
                    watermark='0';
                else
                    watermark='1';
                end
                ExWater(1,counter)=watermark;
            end
            wateronedim(level,:)=ExWater;
            
        end
        %水印恢复
        %extractlevelwatermark=watermarkrestore(wateronedim);
        for level=1:3
            extractlevelwatermark(:,:,level)=gainlevelimage(wateronedim(level,:));
            extractlevelwatermark(:,:,level)=Affine(extractlevelwatermark(:,:,level),a,b,c,d,key,1);
        end
        %toc;
        imwrite(uint8(extractlevelwatermark),'TranslatingAttack extrwater.bmp')
        ncval(k,m)=nc(uint8(extractlevelwatermark),uint8(Water));
        figure(3),subplot(348),imshow(uint8(watermarked)),title('TranslatingAttackcorrected image');
        figure(4),subplot(224),imshow( uint8(extractlevelwatermark)),title(['TranslatingAttack Extracted Watermark nc=',num2str(ncval(k,m))]);
        %------------------------------------------------------
        k=k+1;
    end
end