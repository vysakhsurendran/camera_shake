clc;
clear;
close all;
addpath('NLFMT')

tic

   [fn,pn]=uigetfile('*','Select video');
   obj = VideoReader([pn,fn]);

numFrames = get(obj, 'NumberOfFrames');
  
im1= read(obj,1);
siz=0.3;
im1=imresize(im1,siz);

[row,col,~]=size(im1);
figure,imshow(im1)
pause(0.1)

W=zeros(row,col);
Up=zeros(row,col,3);

for pic=2:60
    pic
    
    im2= read(obj,pic);
    im2=imresize(im2,siz);
    
    for iter=1:3
        Io=im1(:,:,iter);
        Id=im2(:,:,iter);
        
        ptsO  =detectSURFFeatures(Io);
        figure(2),imshow(Io); title('Surf points 1st image')
        hold on;
        plot(ptsO); 
        pause(.5)
        
        ptsD = detectSURFFeatures(Id);
        figure(3),imshow(Id); title('Surf points 2nd image')
        hold on;
        plot(ptsD); 
        pause(.5)
        
        [feaO,validPtsO]  = extractFeatures(Io,ptsO);
        [feaD,validPtsD] = extractFeatures(Id,ptsD);
        
        indP = matchFeatures(feaO,feaD);
        matchedO  = validPtsO(indP(:,1));
        matchedD = validPtsD(indP(:,2));
        
        tform = estimateGeometricTransform(matchedD,matchedO,'similarity');
        
        outputView = imref2d(size(Io));
        Irec  = imwarp(Id,tform,'OutputView',outputView);
        
        
        Im(:,:,iter)=(double(Io)+double(Irec))/2;
    end

% figure,
% subplot(131),imshow(im1);
% subplot(132),imshow(im2);
% subplot(133),imshow(uint8(Im));


% fft

Vi=fft(Im(:));

Vi=reshape(Vi,row,col,3);

W2=sum(abs(Vi),3)/3;

sig=row*col/50;
h=fspecial('gaussian',[5,5],sig);

W1=imfilter(W2,h);

% update weight
W=W+W1;


% update value
Up(:,:,1)=Up(:,:,1)+Vi(:,:,1).*W1;
Up(:,:,2)=Up(:,:,2)+Vi(:,:,2).*W1;
Up(:,:,3)=Up(:,:,3)+Vi(:,:,3).*W1;
% figure,imshow(Up)

end

Up(:,:,1)=Up(:,:,1)./W;
Up(:,:,2)=Up(:,:,2)./W;
Up(:,:,3)=Up(:,:,3)./W;


Up1=abs(ifft(Up(:)));
Up1=reshape(Up1,row,col,3);



% figure,imshow(uint8(Up1))
% pause(0.1)
%====================================================================

%%% NL-means Filter Parameters.
ksize=7;    %%% Neighbor Window Size (should be odd).7
ssize=21;   %%% Search Window Size (should be odd).21
sigmas=5;   %%% Sigma for Gaussian Kernel Generation.5

%%% Wavelet Transform Parameters.
Nlevels=3;
NoOfBands=3*Nlevels+1;
wname='db16'; %% db8 sym8 db16 coif5 bior6.8
sorh='s'; % s or h or t -> trimmed

for iter=1:3
    iter
    xn=Up1(:,:,iter);
    [M,N]=size(xn);
    %%% Noise Level Estimation using Robust Median Estimator.
    if(isequal(wname,'DCHWT'))
        dchw=dchwtf2(xn,1);
        tt1=dchw{1}(:)';
    else
        [ca,ch,cv,cd]=dwt2(xn,wname);
        tt1=cd(:)';
    end
    median_hh2=median(abs(tt1)); %% HH1->Subband containing finest level diagonal details.
    std_dev2=(median_hh2/0.6745);
    
    %%% NL-means Filtering.
    
    im_nl(:,:,iter)=nlmeans_filt2D(xn,sigmas,ksize,ssize,std_dev2);
end

% figure,imshow(uint8(im_nl));

Updenoise=im_nl;

sig=row*col/50;
h=fspecial('gaussian',[3,3],sig);


UpG(:,:,1)=imfilter(Updenoise(:,:,1),h);
UpG(:,:,2)=imfilter(Updenoise(:,:,2),h);
UpG(:,:,3)=imfilter(Updenoise(:,:,3),h);


UpS=2*Updenoise-UpG;


del=0.2;

Up=UpS+del*(Up1-Updenoise);

% figure,
% subplot(131),imshow(uint8(Updenoise));
% subplot(132),imshow(uint8(UpS));
% subplot(133),imshow(uint8(Up));

figure,imshow(uint8(Up));

figure,
subplot(121),imshow(im2);
subplot(122),imshow(uint8(UpS));
disp(abs(psnr(double(im2),UpS)))
disp(abs(mse(double(im2),UpS)))

toc


