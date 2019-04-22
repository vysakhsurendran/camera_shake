clc;
clear;
close all;
addpath('parking_night')
%=========================================================================
im1=imread('REG_001.jpg');
% im1 = imread('IMG_20190227_100104_1.jpg');
siz=0.2;
im1=imresize(im1,siz);

[row,col,~]=size(im1);
figure,imshow(im1)
pause(0.1)

W=zeros(row,col);
Up=zeros(row,col,3);

for pic=2:10
    pic
    if pic<10
        im2=imread(['REG_00' num2str(pic) '.jpg']);
    else
        im2=imread(['REG_0' num2str(pic) '.jpg']);
    end
    
    im2=imresize(im2,siz);
    
    for iter=1:3
        Io=im1(:,:,iter);
        Id=im2(:,:,iter);
        
        
        ptsO  = detectSURFFeatures(Io);
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
        [feaD,validPtsD]  = extractFeatures(Id,ptsD);
        
        indP = matchFeatures(feaO,feaD);
        matchedO  = validPtsO(indP(:,1));
        matchedD = validPtsD(indP(:,2));
        
        tform = estimateGeometricTransform(matchedD,matchedO,'similarity');
        
        outputView = imref2d(size(Io));
        Irec  = imwarp(Id,tform,'OutputView',outputView);
        Im(:,:,iter)=(double(Io)+double(Irec))/2;
    end
        
end
%figure , imshow(Irec);title('Recovered image');