
addpath('bookshelf')

origIm = imread('REG_001.jpg');
littleBlurredIm = imread('REG_005.jpg');
muchBlurredIm = imread('REG_009.jpg');
grayImage1 = rgb2gray(origIm);
grayImage2= rgb2gray(littleBlurredIm);
grayImage3= rgb2gray(muchBlurredIm);

[lpOrigIm,~] =  imgradient(grayImage1);
[lpLittleBlurredIm,~] = imgradient(grayImage2);
[lpMuchBlurredIm,~] = imgradient(grayImage3);



% Number of pixels to look at: 0.1%
nPx = round(0.001*numel(origIm));

% Sort values to pick top values
sortedOrigIm = sort(lpOrigIm(:));
sortedLittleBlurredIm = sort(lpLittleBlurredIm(:));
sortedMuchBlurredIm = sort(lpMuchBlurredIm(:));

% Calculate measure
measureOrigIm = median(sortedOrigIm(end-nPx+1:end));
measureLittleBlurredIm = median(sortedLittleBlurredIm(end-nPx+1:end));
measureMuchBlurredIm = median(sortedMuchBlurredIm(end-nPx+1:end));