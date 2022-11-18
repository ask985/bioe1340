% assignment 3 part 1
% 10/25/22
clear; clc; clear all;

%% part 1 question 1a
% average the snr, 2nd is averaging the images and then getting the snr
imgs = dir("RatBrain/*.TIFF");
imgs = {imgs.name}';
img = zeros(1200,1600,numel(imgs));
snr = [];

I2 = imread('RatBrain/I_VOL_1.TIFF');
imshow(I2);
%drawfreehand to select yellow region by hand
fh = drawfreehand;
mask = createMask(fh);

%% 1a continued
for i=1:numel(imgs)
    img(:,:,i) = im2double(rgb2gray(imread(['RatBrain/', imgs{i}]))); %1200x1600x3
    %crop out everything except infarction (yellow pixel) region using the
    %mask i made previously
    roi = img(:,:,i);
    roi = roi(mask==1);

    sgnl = mean(roi(roi>0));
    noise = std(roi(roi>0));
    %calculate the snr of each image
    snr(i) = 10 * log10(sgnl/noise);
end

%take the average of the SNRs we calculated for each image
avesnr = sum(snr)/numel(snr); %average snr of the 50 images is 8.02


%% 1b

img_avg = mean(img,3);
roi_avg = img_avg;
roi_avg = roi_avg(mask==1);

sgnl = mean(roi_avg(roi_avg>0));
noise = std(roi_avg(roi_avg>0));
snr_avg = 10 * log10(sgnl/noise); %9.45 now - it's 17.8% higher 
%more smooth because we lowered the standard deviation so snr can be higher

%for loop. see relationship between sample size and snr
snr_box = [];
for i=1:50
    img_temp = img(:,:,1:i);
    avg_temp = mean(img_temp,3);
    roi_avg2 = avg_temp;
    roi_avg2 = roi_avg2(mask==1);
    
    sgnl2 = mean(roi_avg2(roi_avg2>0));
    noise2 = std(roi_avg2(roi_avg2>0));
    
    snr_box(i) = 10 * log10(sgnl2/noise2);
    
end
figure;
plot(sqrt(1:50), snr_box); %this gives a square root plot, so the 
xlabel('sqrt(N)');
ylabel('SNR');


%% part 1 part 2 (CNR) creating yellow mask
I2 = imread('RatBrain/I_VOL_1.TIFF');
imshow(I2);
fh1 = drawfreehand;
ymask = createMask(fh1);

%% making red mask
fh2 = drawfreehand;
rmask = createMask(fh2);

%% part 1 part 2 continued
imgs = dir("RatBrain/*.TIFF");
imgs = {imgs.name}';
%img = zeros(1200,1600,numel(imgs));

slice = [];
cnrs = [];
for i = 1 : numel(imgs)
    slice(:,:,i) = im2double(rgb2gray((imread(['RatBrain/', imgs{i}]))));
    
    yroi = slice(:,:,i);
    yroi = yroi(ymask==1);
    rroi = slice(:,:,i);
    rroi = rroi(rmask==1);

    %find cnr
    ysgnl = mean(yroi(yroi>0));
    ynoise = std(yroi(yroi>0));
    rsgnl = mean(rroi(rroi>0));
    rnoise = std(rroi(yroi>0));

    cnr = (ysgnl - rsgnl)/sqrt((ynoise^2 + rnoise^2)/2);
    cnrs(i) = cnr;

end

%average cnr of individual images
cnr1 = sum(cnrs)/numel(cnrs) %average cnr of indiv images: 2.4367

%% part 1 part 2 continued 
slice_avg = mean(slice,3);

yroi_avg = slice_avg; %2d slice
yroi_avg = yroi_avg(ymask==1);
rroi_avg = slice_avg;
rroi_avg = rroi_avg(rmask==1);

%find cnr
ysgnl_avg = mean(yroi_avg(yroi_avg>0));
ynoise_avg = std(yroi_avg(yroi_avg>0));
rsgnl_avg = mean(rroi_avg(rroi_avg>0));
rnoise_avg = std(rroi_avg(yroi_avg>0));

cnr_avg = (ysgnl_avg - rsgnl_avg)/sqrt((ynoise_avg^2 + rnoise_avg^2)/2)
%cnr_avg = 3.3722, which is an 38.39% increase.


%% part 1 question 2
%reading shepp logan phantom
I = phantom('Modified Shepp-Logan');
imshow(I);

%loads variable kspaceSignal
load('2DSheppLoganKSpace.mat');

%take inverse fourier transform of kspaceSignal
data = ifft2(kspaceSignal);
I2 = fftshift(data);
I2_comb = abs(I2);
minI2 = min(I2_comb, [], "all");
maxI2 = max(I2_comb, [], "all");

%% using imshow
%figure 1
imshow(I2_comb, [minI2, maxI2]);
%set(gcf,'position',[50,50,600,600])

%% using imagesc
figure;
%figure 2
imagesc(I2_comb);
colormap("gray")
set(gcf,'position',[50,50,600,600])
%more flexible - will automatically adjust image dimensions

%% visualize k-space, centering
figure;
%center it
%flips the image from center to bounday
Iq3 = ifftshift(kspaceSignal); 

Iq3 = abs(Iq3);
%figure 3
imagesc(Iq3); %visualizes k space
%original center now located at the four corners of the image

%% for part2 question 1
%real and imaginary parts of I2
I2_real = real(I2);
I2_imag = imag(I2);

figure; %plot of real
%figure 4
imagesc(I2_real);

figure; %plot of imaginary
%figure 5
imagesc(I2_imag);

%real and imaginary components are orthogonal

%% original k-space visualization
figure;
Iq4 = kspaceSignal;
Iq4 = abs(Iq4);
imagesc(Iq4);
%this figure shows that the original k-space is already centered