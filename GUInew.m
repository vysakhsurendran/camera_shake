function varargout = GUInew(varargin)
% GUINEW MATLAB code for GUInew.fig
%      GUINEW, by itself, creates a new GUINEW or raises the existing
%      singleton*.
%
%      H = GUINEW returns the handle to a new GUINEW or the handle to
%      the existing singleton*.
%
%      GUINEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUINEW.M with the given input arguments.
%
%      GUINEW('Property','Value',...) creates a new GUINEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUInew_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUInew_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUInew

% Last Modified by GUIDE v2.5 25-Apr-2019 19:39:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUInew_OpeningFcn, ...
                   'gui_OutputFcn',  @GUInew_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUInew is made visible.
function GUInew_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUInew (see VARARGIN)

% Choose default command line output for GUInew
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUInew wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUInew_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.*';'*.bmp';'*.jpg';'*.gif'}, 'Select the reference image');
imm = imread([pathname,filename]);
siz=0.4;
imm=imresize(imm,siz);
[row,col,~]=size(imm);
axes(handles.axes1);
imshow(imm);
%handles.pathname1 =ab
original=imread([pathname,filename]);

%blur level estimation of input image first method
%{
I = double(original);
[y x] = size(I);
Hv = [1 1 1 1 1 1 1 1 1]/9;
Hh = Hv';
B_Ver = imfilter(I,Hv);%blur the input image in vertical direction
B_Hor = imfilter(I,Hh);%blur the input image in horizontal direction
D_F_Ver = abs(I(:,1:x-1) - I(:,2:x));%variation of the input image (vertical direction)
D_F_Hor = abs(I(1:y-1,:) - I(2:y,:));%variation of the input image (horizontal direction)
D_B_Ver = abs(B_Ver(:,1:x-1)-B_Ver(:,2:x));%variation of the blured image (vertical direction)
D_B_Hor = abs(B_Hor(1:y-1,:)-B_Hor(2:y,:));%variation of the blured image (horizontal direction)
T_Ver = D_F_Ver - D_B_Ver;%difference between two vertical variations of 2 image (input and blured)
T_Hor = D_F_Hor - D_B_Hor;%difference between two horizontal variations of 2 image (input and blured)
V_Ver = max(0,T_Ver);
V_Hor = max(0,T_Hor);
S_D_Ver = sum(sum(D_F_Ver(2:y-1,2:x-1)));
S_D_Hor = sum(sum(D_F_Hor(2:y-1,2:x-1)));
S_V_Ver = sum(sum(V_Ver(2:y-1,2:x-1)));
S_V_Hor = sum(sum(V_Hor(2:y-1,2:x-1)));
blur_F_Ver = (S_D_Ver-S_V_Ver)/S_D_Ver;
blur_F_Hor = (S_D_Hor-S_V_Hor)/S_D_Hor;
blur = max(blur_F_Ver,blur_F_Hor);
set(handles.text5, 'string',blur);
%}

%blur level estimation of input image second method
%{
origIm = imread([pathname,filename]);
grayImage1 = rgb2gray(origIm);
[lpOrigIm,~] =  imgradient(grayImage1);
% Number of pixels to look at: 0.1%
nPx = round(0.001*numel(origIm));
% Sort values to pick top values
sortedOrigIm = sort(lpOrigIm(:));
% Calculate measure
c = median(sortedOrigIm(end-nPx+1:end));
%}
handles.fi1 = filename;
handles.pa1 = pathname;
guidata(hObject,handles);
%set(handles.text5, 'string', c);
% set(handles.edit3,'string',Mean);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%a3=handles.pathname1;
I1 = handles.fi1;
I2 = handles.pa1;


addpath(I2)
im1=imread(I1);
%im1 = imread('IMG_20190227_100104_1.jpg');
siz=0.4;
im1=imresize(im1,siz);

[row,col,~]=size(im1);
%figure,imshow(im1)
%pause(0.1)

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
        %figure(2),imshow(Io); title('Surf points 1st image')
        %hold on;
        %plot(ptsO); 
        %pause(.5)

        ptsD = detectSURFFeatures(Id);
        %figure(3),imshow(Id); title('Surf points 2nd image')
        %hold on;
        %plot(ptsD); 
        %pause(.5)
        
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



%figure,imshow(uint8(Up1))
%pause(0.1)
%====================================================================

%%% NL-means Filter Parameters.
ksize=7;    %%% Neighbor Window Size (should be odd).7
ssize=21;   %%% Search Window Size (should be odd).21
sigmas=5;   %%% Sigma for Gaussian Kernel Generation.5

%%% Wavelet Transform Parameters.
Nlevels=3;
NoOfBands=3*Nlevels+1;
wname='db8'; %% db8 sym8 db16 coif5 bior6.8
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

%figure,imshow(uint8(im_nl));

Updenoise=im_nl;

sig=row*col/50;
h=fspecial('gaussian',[3,3],sig);


UpG(:,:,1)=imfilter(Updenoise(:,:,1),h);
UpG(:,:,2)=imfilter(Updenoise(:,:,2),h);
UpG(:,:,3)=imfilter(Updenoise(:,:,3),h);


UpS=2*Updenoise-UpG;


del=0.2;

Up=UpS+del*(Up1-Updenoise);

%figure,
%subplot(131),imshow(uint8(Updenoise));
%subplot(132),imshow(uint8(UpS));
%subplot(133),imshow(uint8(Up));

%figure,imshow(uint8(Up));

%figure,
%subplot(121),imshow(im1);
%subplot(122),imshow(uint8(UpS));
axes(handles.axes3);
imshow(uint8(UpS));
%{
origIm = imread(uint8(UpS));
grayImage1 = rgb2gray(uint8(UpS);
[lpOrigIm,~] =  imgradient(grayImage1);
% Number of pixels to look at: 0.1%
nPx = round(0.001*numel(uint8(UpS)));
% Sort values to pick top values
sortedOrigIm = sort(lpOrigIm(:));
% Calculate measure
cb = median(sortedOrigIm(end-nPx+1:end));
%textLabel = fprintf('%f', cb);
set(handles.text6, 'string', cb);
%}

imwrite(uint8(UpS),fullfile(I2,'output.jpg'))
%psnr_=getPSNR(uint8(im1),uint8(UpS));
disp(abs(psnr(double(im2),UpS)))
ssim_=getMSSIM(uint8(UpS),uint8(im1));
%fprintf('PSNR= %f - SSIM= %f\n',psnr_,ssim_);
fprintf('SSIM=%f\n',ssim_);
set(handles.text6, 'string', ssim_);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)`

close all;
