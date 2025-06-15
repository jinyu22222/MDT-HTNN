

%e.g., img_house or mri
X0=data;

missing_rate=0.4; % sample ratio(SR)=1-missing_rate
tau=[10,10,1]; % delay window size
[height, width, band] = size(X0);             
dim = [height, width, band];
miss_num=prod(dim)*missing_rate;
sampling_rate = 1-missing_rate;
m          = round(prod(dim)*sampling_rate);
sort_dim   = randperm(prod(dim));
Omega      = sort_dim(1:m); % sampling pixels' index
Omega_MDT_HTNN= zeros(dim);
Omega_MDT_HTNN(Omega) = 1; % observed Img
 methodName='MDT-HTNN';

tic
[X, err,iter]=MDT_TNN_DCT(X0, tau, Omega_MDT_HTNN, miss_num);
Time=toc;

%% Show result
fprintf('\n');    
fprintf('================== QA Results =====================\n');
fprintf(' %8.8s    %5.5s    %5.5s    %5.5s     %5.5s  \n',...
    'Method', 'MPSNR', 'MSSIM', 'MFSIM',  'Time');
    fprintf(' %8.8s   %5.3f    %5.3f    %5.3f    %5.3f   \n',...
        methodName, PSNR, SSIM, FSIM, Time);
fprintf('================== Show Results =====================\n');

imshow(X)


