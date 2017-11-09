%% Question 1.1 %%
clear all; close all; clc;

I = imread('..\Tiffany.jpg');
figure(1);
subplot(1,2,2); imshow(I,[]);
title('Original Picture','fontsize',16);
I_fourier_tiffany = fftshift(fft2(I));
subplot(1,2,1); imshow(log10(abs(I_fourier_tiffany)+1),[]);
title('Original Picture in Frequency Domain','fontsize',16);

%% Question 1.2 %%
I_LPF = I_fourier_tiffany;
I_LPF([1:242,269:512],[1:242,269:512]) = 0;
figure(2);
subplot(1,2,1); imshow(log10(abs(I_LPF)+1),[]);
title('Picture after LPF in Frequency Domain','fontsize',16);
I_reversed = ifft2(ifftshift(I_LPF));
subplot(1,2,2); imshow(abs(I_reversed),[]);
title('Picture after LPF','fontsize',16);

%% Question 1.6 %%
I_fourier_sum1 = sum(abs(I_fourier_tiffany));
[~,vec1] = sort(I_fourier_sum1);
vec1 = vec1((512-26+1):512);
I_fourier_sum2 = sum(abs(I_fourier_tiffany'));
[~,vec2] = sort(I_fourier_sum2);
vec2 = vec2((512-26+1):512);

I_most_energy_rows_columns = I_fourier_tiffany;
I_most_energy_rows_columns(setdiff(1:512,vec2),setdiff(1:512,vec1)) = 0; % sort according to rows and columns

figure(3);
subplot(1,2,1); imshow(log10(abs(I_most_energy_rows_columns)+1),[]);
title('Picture after most energy taken in rows and columns Frequency Domain','fontsize',12);
I_reversed = ifft2(ifftshift(I_most_energy_rows_columns));
subplot(1,2,2); imshow(abs(I_reversed),[]);
title('Picture after most energy taken in rows and columns','fontsize',12);

%% Question 1.7 %%
I_most_energy = I_fourier_tiffany;
I_most_energy_abs = abs(I_most_energy);
matrix_cells_sorted = sort(I_most_energy_abs(:));
threshold_val = matrix_cells_sorted(512*512-13107);
I_most_energy(I_most_energy_abs<threshold_val) = 0; % zero not max energy pixels

figure(4);
subplot(1,2,1); imshow(log10(abs(I_most_energy)+1),[]);
title('Picture after most energy taken in Frequency Domain','fontsize',16);
I_reversed = ifft2(ifftshift(I_most_energy));
subplot(1,2,2); imshow(abs(I_reversed),[]);
title('Picture after most energy taken','fontsize',16);

%% Question 2.1 %%
clear all; close all; clc;
I_brad = imread('..\brad.jpg');
I_baboon = imread('..\baboon.jpg');

figure(1);
subplot(1,2,1); imshow(I_brad,[]);
title('Original Picture','fontsize',16);
I_fourier_brad = fftshift(fft2(I_brad));
I_fourier_brad_amp = abs(I_fourier_brad);
I_fourier_brad_phase = angle(I_fourier_brad);
subplot(1,2,2); imshow(log10(abs(I_fourier_brad)+1),[]);
title('Original Picture in Frequency Domain','fontsize',16);

figure(2);
subplot(1,2,1); imshow(I_baboon,[]);
title('Original Picture','fontsize',16);
I_fourier_baboon = fftshift(fft2(I_baboon));
I_fourier_baboon_amp = abs(I_fourier_baboon);
I_fourier_baboon_phase = angle(I_fourier_baboon);
subplot(1,2,2); imshow(log10(abs(I_fourier_baboon)+1),[]);
title('Original Picture in Frequency Domain','fontsize',16);

%% Question 2.2 %%
% crossing amplitudes and phase of 2 images
I_brad_amp_baboon_phase = I_fourier_brad_amp.*exp(1i.*I_fourier_baboon_phase);
I_brad_phase_baboon_amp = I_fourier_baboon_amp.*exp(1i.*I_fourier_brad_phase);

figure(3);
subplot(1,2,1); imshow(abs(ifft2(ifftshift(I_brad_amp_baboon_phase))),[]);
title('Brad Pitt Amplitude, Baboon Phase','fontsize',16);
subplot(1,2,2); imshow(abs(ifft2(ifftshift(I_brad_phase_baboon_amp))),[]);
title('Brad Pitt Phase, Baboon Amplitude','fontsize',16);

%% Question 2.3 %%
% creating images with random phase or amplitude
I_brad_amp_random_phase = I_fourier_brad_amp.*exp(1i.*(rand(512)-0.5)*2*pi);
I_brad_phase_random_amp = randi([0,255],512).*exp(1i*I_fourier_brad_phase);

figure(4);
subplot(1,2,1); imshow(abs(ifft2(ifftshift(I_brad_amp_random_phase))),[]);
title('Brad Pitt Amplitude, Random Phase','fontsize',16);
subplot(1,2,2); imshow(abs(ifft2(ifftshift(I_brad_phase_random_amp))),[]);
title('Brad Pitt Phase, Random Amplitude','fontsize',16);

%% Question 3.1 %%
clear all; close all; clc;
r_norm = linspace(0,1,256);
r = r_norm*255;
gamma = [0.04, 0.1, 0.2, 0.4, 0.67, 1, 1.5, 2.5, 5, 10, 25]; 

s_norm = zeros(length(gamma),length(r_norm));

figure(1);
myMap = rand(length(gamma), 3);

for i=1:length(gamma)
    s_norm(i,:) = r_norm.^gamma(i); % compute gamma correction values
    plot(r, 255*s_norm(i,:),'Color',myMap(i,:)); hold on;
end

title('s(r) function for different gamma values','fontsize',16);
legend('\gamma = 0.04', '\gamma = 0.1', '\gamma = 0.2', '\gamma = 0.4', '\gamma = 0.67', ...
    '\gamma = 1', '\gamma = 1.5', '\gamma = 2.5', '\gamma = 5', '\gamma = 10', '\gamma = 25') 
ylabel('s(r)','fontsize',16); xlabel('r','fontsize',16);

%% Question 3.2 %%
I = imread('..\mri_spine.jpg');
LUT_1 = uint8(round(255*s_norm(6,:)/max(s_norm(6,:))));
LUT_2 = uint8(round(255*s_norm(5,:)/max(s_norm(5,:))));
figure(2);
hist_org = imhist(I);
subplot(2,2,1); imshow(I,[]); title('Original Picture','fontsize',16);
subplot(2,2,3); imhist(LUT_1(double(I)+1)); title('Histogram of Original Picture','fontsize',16);

s_norm = im2uint8(s_norm);

% create gamma corrected image
I_corrected = intlut(I,s_norm(5,:));
subplot(2,2,2); imshow(I_corrected,[]); title('Corrected Picture with \gamma = 0.67','fontsize',16);

histogram_corrected = imhist(I_corrected);
subplot(2,2,4); imhist(LUT_2(double(I)+1)); title('Histogram of Corrected Picture with \gamma = 0.67','fontsize',16);

%% Question 3.3 %%
% Full histogram equalization
[I_equal, LUT_equal] = histogram_equal(I,0);
figure(3);
subplot(1,3,1); imhist(I_equal); title('Histogram of Equalized Image','fontsize',16);
subplot(1,3,2); stairs(LUT_equal); title('LUT of Equalized Image','fontsize',16);
subplot(1,3,3); imshow(I_equal); title('Equalized Image','fontsize',16);

%% Question 3.4 %%
% Histogram Equalization above a threshold
[I_equal_2, LUT_equal_2] = histogram_equal(I,3);
figure(4);
subplot(1,3,1); imhist(I_equal_2); title('Histogram of Equalized Image - With threshold=3','fontsize',12);
subplot(1,3,2); stairs(LUT_equal_2); title('LUT of Equalized Image - With threshold=3','fontsize',12);
subplot(1,3,3); imshow(I_equal_2); title('Equalized Image - With threshold=3','fontsize',12);

%% Question 3.5 %%
% dividing picture to 4 areas according to intensity distribution
I1 = I(:,1:136);
I2 = I(:,137:216);
I3 = I(:,217:295);
I4 = I(:,296:373);

[I_equal_I1, LUT_equal_I1] = histogram_equal(I1,8);
[I_equal_I2, LUT_equal_I2] = histogram_equal(I2,0);
[I_equal_I3, LUT_equal_I3] = histogram_equal(I3,10);
[I_equal_I4, LUT_equal_I4] = histogram_equal(I4,6);

% re-attach all 4 parts of the image
I_equal_final = [I_equal_I1, I_equal_I2, I_equal_I3, I_equal_I4];
figure(5);
subplot(1,4,1); imhist(I1); title('Histogram of Original Partial Image #1','fontsize',8);
subplot(1,4,2); imhist(I2); title('Histogram of Original Partial Image #2','fontsize',8);
subplot(1,4,3); imhist(I3); title('Histogram of Original Partial Image #3','fontsize',8);
subplot(1,4,4); imhist(I4); title('Histogram of Original Partial Image #4','fontsize',8);

figure(6);
subplot(1,4,1); imhist(I_equal_I1); title('Histogram of Equalized Image #1 - With threshold=8','fontsize',8);
subplot(1,4,2); imhist(I_equal_I2); title('Histogram of Equalized Image #2 - With threshold=0','fontsize',8);
subplot(1,4,3); imhist(I_equal_I3); title('Histogram of Equalized Image #3 - With threshold=10','fontsize',8);
subplot(1,4,4); imhist(I_equal_I4); title('Histogram of Equalized Image #4 - With threshold=6','fontsize',8);

figure(7);
imshow(I_equal_final,[]);
title('Final Processed Image - after merging 4 parts','fontsize',12);