
%
% Demo TV/L1 solve
%

clc,clear; % close all; 
% path(path,genpath(pwd));

% path(path,'./WNNM_code/');
path(path,'./Noisy_data/');
path(path,'./Blurred_noisy_data/');

% addpath('E:\impulsive-noise\code\');

%% load original image


 
 I = double(imread('parrot.png'))/255;

 % I = double(imread('house256.png'))/255;
 
 % I = double(imread('bridge.tiff'))/255;

% I = double(rgb2gray(imread('building_org.png')))/255;










 H = fspecial('average',1);

% H = fspecial('gaussian',7,5);

% level = 0.7;
% Bn = imfilter(I,H,'circular');
% Bn = imnoise(Bn,'salt & pepper',level);







 tol = 1e-6;
 iter = 4000;
 rho = 14;


%% regularization parameter
 
mu = [1,9,10,7,6]; % for parrot denoising

% mu = [10,9,15,13,3]; % for bridge denoising

% mu = [20,15,7,6,1]; % for house denoising

% mu = [100,110,140,180,200]; % for building denoising

%  mu = [1,2,4,5,12]; % for parrot deblurring

% mu = [4,5,5,5,5]; % for house deblurring 

% mu = [2,5,7,12,18]; % for bridge deblurring 

% mu = [6,8,14,26,50]; % for building deblurring




%%


% [x,PSNR,mssim,i] = PD_Constrained_TVL1_three(I,fb,H,tau,sigma,tol,iter,N,mu);
% [x,PSNR,mssim,i] = PD_Constrained_TVL1_two(I,Bn,fb,H,tau,sigma,tol,iter,N,mu,mu1);

% [ x,k,SNR,PSNR,SIM,tg] = ADMM_impulsive(I,fb,H,iter,tol,mu,N);
% [ x,k,SNR,PSNR,SIM,tg] = Variant_ADMM_impulsive(I,fb,H,iter,tol,mu,N);


psnr_input = zeros(1,length(mu));
ssim_input = zeros(1,length(mu));
x1 = cell(1,length(mu));
x2 = cell(1,length(mu));
x3 = cell(1,length(mu));
k1 = zeros(1,length(mu));
k2 = zeros(1,length(mu));
k3 = zeros(1,length(mu));
SNR1 = zeros(1,length(mu));
SNR2 = zeros(1,length(mu));
SNR3 = zeros(1,length(mu));
PSNR1 = zeros(1,length(mu));
PSNR2 = zeros(1,length(mu));
PSNR3 = zeros(1,length(mu));
SIM1 = zeros(1,length(mu));
SIM2 = zeros(1,length(mu));
SIM3 = zeros(1,length(mu));

tg1 = cell(1,length(mu));
tg2 = cell(1,length(mu));
tg3 = cell(1,length(mu));

PSNR1k = cell(1,length(mu));
PSNR2k = cell(1,length(mu));
PSNR3k = cell(1,length(mu));


for j = 1:length(mu)
    
% storageName = strcat('parrot',num2str(j),'.mat');

 noise_image = ['parrot',num2str(j),'.mat'];
% noise_image = ['house',num2str(j),'.mat'];
% noise_image = ['bridge',num2str(j),'.mat'];
% noise_image = ['building',num2str(j),'.mat'];

% noise_image = ['parrot_GS_7_5',num2str(j),'.mat'];
% noise_image = ['house_GS_7_5',num2str(j),'.mat'];
% noise_image = ['bridge_GS_7_5',num2str(j),'.mat'];
% noise_image = ['building_GS_7_5',num2str(j),'.mat'];


load(noise_image);

psnr_input(j) = psnr(Bn,I);
ssim_input(j) = ssim(Bn,I);


%% detect the location of noisy images

[m,n] = size(I); 
N = [];

for i1 = 1:m
    for j1 = 1:n
        if Bn(i1,j1) == 0 || Bn(i1,j1) == 1
            N(i1,j1) = 0;
        else
            N(i1,j1) = 1;
        end
    end
end

fb = Bn.*N;  % input image
 


[ x1{j},k1(j),SNR1(j),PSNR1(j),SIM1(j),tg1{j},PSNR1k{j}] = Dong_ADMM_constrained_impulsive(I,fb,H,iter,tol,N,rho);
 

% imwrite(x1{j},['F:\Impulse_noise\Restored_images_denoising\Dong_parrot',num2str(j),'.png']);
% imwrite(x1{j},['F:\Impulse_noise\Restored_images_denoising\Dong_house',num2str(j),'.png']);
% imwrite(x1{j},['F:\Impulse_noise\Restored_images_denoising\Dong_bridge',num2str(j),'.png']);
% imwrite(x1{j},['F:\Impulse_noise\Restored_images_denoising\Dong_building',num2str(j),'.png']);

% imwrite(x1{j},['F:\Impulse_noise\Restored_images_deblurring\Dong_deblurring_parrot',num2str(j),'.png']);
% imwrite(x1{j},['F:\Impulse_noise\Restored_images_deblurring\Dong_deblurring_house',num2str(j),'.png']);  
% imwrite(x1{j},['F:\Impulse_noise\Restored_images_deblurring\Dong_deblurring_bridge',num2str(j),'.png']); 
% imwrite(x1{j},['F:\Impulse_noise\Restored_images_deblurring\Dong_deblurring_building2',num2str(j),'.png']); 

  [ x2{j},k2(j),SNR2(j),PSNR2(j),SIM2(j),tg2{j},PSNR2k{j}] = Box_Dong_ADMM_constrained_impulsive(I,fb,H,iter,tol,N,rho);
  
% restored_image2 = ['Constr_Dong_parrot',num2str(j)];
% imwrite(x2{j},['F:\Impulse_noise\Restored_images_denoising\Constr_Dong_parrot',num2str(j),'.png']);
% imwrite(x2{j},['F:\Impulse_noise\Restored_images_denoising\Constr_Dong_house',num2str(j),'.png']);
%  imwrite(x2{j},['F:\Impulse_noise\Restored_images_denoising\Constr_Dong_bridge',num2str(j),'.png']);
% imwrite(x2{j},['F:\Impulse_noise\Restored_images_denoising\Constr_Dong_building',num2str(j),'.png']);

% imwrite(x2{j},['F:\Impulse_noise\Restored_images_deblurring\Constr_Dong_deblurring_parrot',num2str(j),'.png']);
% imwrite(x2{j},['F:\Impulse_noise\Restored_images_deblurring\Constr_Dong_deblurring_house',num2str(j),'.png']);
% imwrite(x2{j},['F:\Impulse_noise\Restored_images_deblurring\Constr_Dong_deblurring_bridge',num2str(j),'.png']);
% imwrite(x2{j},['F:\Impulse_noise\Restored_images_deblurring\Constr_Dong_deblurring_building2',num2str(j),'.png']);

 [ x3{j},k3(j),SNR3(j),PSNR3(j),SIM3(j),tg3{j},PSNR3k{j}] = ADMM_constrained_impulsive(I,fb,H,iter,tol,mu(j),N,rho);
 
% restored_image3 = ['Ours_parrot',num2str(j)];
% imwrite(x3{j},['F:\Impulse_noise\Restored_images_denoising\Ours_parrot',num2str(j),'.png']);
% imwrite(x3{j},['F:\Impulse_noise\Restored_images_denoising\Ours_house',num2str(j),'.png']);
% imwrite(x3{j},['F:\Impulse_noise\Restored_images_denoising\Ours_bridge',num2str(j),'.png']);
% imwrite(x3{j},['F:\Impulse_noise\Restored_images_denoising\Ours_building',num2str(j),'.png']);

% imwrite(x3{j},['F:\Impulse_noise\Restored_images_deblurring\Ours_deblurring_parrot',num2str(j),'.png']);
% imwrite(x3{j},['F:\Impulse_noise\Restored_images_deblurring\Ours_deblurring_house',num2str(j),'.png']); 
% imwrite(x3{j},['F:\Impulse_noise\Restored_images_deblurring\Ours_deblurring_bridge',num2str(j),'.png']);
% imwrite(x3{j},['F:\Impulse_noise\Restored_images_deblurring\Ours_deblurring_building2',num2str(j),'.png']);


end


