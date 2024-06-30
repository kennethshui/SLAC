function [mask, RfacR, RfacF, RFD]= HIO_ER (F2D,supp,iter,beta,showim,modelimg) 
%function [mask, RfacR, RfacF, RFD]= HIO_ER (F2D,supp,iter,beta,showim,modelimg)
%Phase retrieval using HIO+ER, based on 'OSS_revised.m' on 07/03/22.
%F2D:           Diffraction pattern of the oversampled vesicle image, with the central speckle removed, and 5% Poisson noise (Rnoise) added as defined above.
%supp [52, 28]: Size of rectangular support along the first (y) dimension and second (x) dimension, always placed at the center of image.
%iter [20]:     Total number of iterations performed by the HIO+ER. Inside each iteration, 400HIO + 200ER.
%beta [0.9]:    Parameter used in a typical HIO reconstruction (see Fienup, Opt. Let. 1978).
%showim [1]:    Defines whether or not images showing progress are displayed (1 is shown, 0 not shown).
%modelimg:      Image of the oversampled vesicle model used for the OSS reconstruction.
%mask:          An image of the mask used as the support in the OSS reconstruction.
%RfacR:         Error measurement in real space, defined above, computed alongside Rf, 1000x1 vector.
%RfacF:         Error measurement in Fourier space, defined above, computed every 2 iterations, 1000x1 vector.
%RFD:           Final result image obtained from the OSS reconstruction.
%Usage:         load model.mat; load Ten_perc_noise_Diffraction_Pattern.mat;
%               [mask, RfacR, RfacF, RFD]= HIO_ER(diffpat,[60,30],2000,0.9,1,model,0);

close all;

%% General assignment of variables
[Rsize,Csize] = size(F2D)
R2D=zeros(Rsize,Csize,1,'single');    %10 frame of 2-D real sapce data.
realgroup=zeros(Rsize,Csize,2,'single');
realgroup(:,:,1)=modelimg;
HIO_iter = 400*1;
ER_iter = 200;

 
%% Define support
Rcenter = ceil(Rsize/2);
Ccenter = ceil(Csize/2);
Rsupport = supp(1);
Csupport = supp(2);
half_Rsupport = ceil(Rsupport/2);
half_Csupport = ceil(Csupport/2);
support = zeros(Rsize,Csize,'single');
%support(Rcenter-half_Rsupport+1:Rcenter+half_Rsupport-1,Ccenter-half_Csupport+1:Ccenter+half_Csupport-1) = 1;  %Original support code.
support(Rcenter-half_Rsupport+1:Rcenter+half_Rsupport,Ccenter-half_Csupport+1:Ccenter+half_Csupport) = 1;  %Modified support code.
mask=support;

%Inverse support
inv_support = ones(Rsize,Csize,'single');
inv_support(Rcenter-half_Rsupport+1:Rcenter+half_Rsupport,Ccenter-half_Csupport+1:Ccenter+half_Csupport) = 0;

figure(10);
imagesc(mask);
figure(11);
imagesc(inv_support);

%% Draw the real space image with support region.
figure(3)
imagesc(modelimg);
rectangle('Position',[Ccenter-half_Csupport+1, Rcenter-half_Rsupport+1, Csupport, Rsupport], 'EdgeColor', 'r', 'LineWidth', 1);
title('Model input image with support region');
 

 
%% Assign random phases
rng('shuffle','twister');
phase_angle = rand(Rsize,Csize,'single'); 
 
%% Define initial k, r space
initial_k=F2D; initial_k(initial_k==-1)=0;
k_space = initial_k.*exp(1i*phase_angle);
buffer_r_space = single(real(ifftn(ifftshift(k_space))));    
 
%% Preallocate error arrays
RfacF = zeros(ceil(iter/2),1,'single');  counter1=0; 
RfacR = zeros(ceil(iter/2),1,'single');  counter2=0;
 
%% Image display argument
if showim==1    
    figure(1),
end

 
%% OSS iterations
for iteration = 1:iter

    %% HIO with Support & Positivity constraint
    for iteration_HIO = 1 : HIO_iter  
        
        %Real space and support constraint.
        r_space = real(ifftn(ifftshift(k_space)));
        sample = r_space.*mask;    %Support constraint.
        r_space = buffer_r_space-beta*r_space;
        %sample(sample<0)=r_space(sample<0);  %Positivity constraint.
        r2 = r_space.*inv_support;
        r_space = sample + r2;
        
        
        
        % K-space magnitude constraintn
        buffer_r_space = r_space;
        k_space = fftshift(fftn(r_space));  
        phase_angle = angle(k_space);        
        k_space = F2D.*exp(1i*phase_angle);   %Impose K-space magnitude here.
    end
    
    %% ER with Support & Positivity constraint
    for iteration_ER = 1 : ER_iter 
        
        %Real space and support constraint.
        r_space = real(ifftn(ifftshift(k_space)));
        sample = r_space.*mask;    %Support constraint.
        %sample(sample<0)=r_space(sample<0);  %Positivity constraint.
        r_space = sample;
        
        
        % K-space magnitude constraint 
        k_space = fftshift(fftn(r_space));  
        phase_angle = angle(k_space);        
        k_space = F2D.*exp(1i*phase_angle);   %Impose K-space magnitude here.
        
    end
        
                   
        %% Calculate errors, computed every 2 iterations   
        if rem(iteration,2)==0  
            
            %% Calculate error in reciprocal(Fourier) space
            Ktemp = sample;
            Ktemp = abs(fftshift(fftn(Ktemp)));
            errorF = sum(sum(abs(Ktemp(F2D~=-1)-F2D(F2D~=-1)))) / sum(sum(abs(F2D(F2D~=-1))));    %F2D is the input diffraction pattern magnitude.
            counter1=counter1+1; RfacF(counter1) = errorF;
            
            realgroup(:,:,2)=sample;
            realgroup2=realgroup(Rcenter-half_Rsupport-1:Rcenter+half_Rsupport+1,Ccenter-half_Csupport-1:Ccenter+half_Csupport+1,:);
            realgroup2 = align2(realgroup2,0,iteration);
            
            %% Figure shows progress
            if showim==1
                figure (1),
                
                subplot(2,2,1), imagesc(squeeze(realgroup2(:,:,1))), axis image, title('HIO');              
                subplot(2,2,2), imagesc(squeeze(realgroup2(:,:,2))), axis image, title(strcat('Iteration = ', int2str(iteration)));
                subplot(2,2,3), plot(RfacF), axis([0 ceil(iteration/2) 0 0.8]), title(strcat('F error=',int2str(errorF*100),'%')); grid on;
                subplot(2,2,4),  grid on;
                drawnow
            end
        end
end
 

%% Show image: sum of best 4 steps

RFD = r_space;

if showim==1
figure(2),
subplot(2,2,1), imagesc(squeeze(RFD)), axis image, title('Retrieved image at Real Space'); 
subplot(2,2,2), imagesc(squeeze(mask)), axis image, title('Mask/support area');
subplot(2,2,3), plot(RfacF), axis tight, title('F-space error'); grid on;
subplot(2,2,4), grid on;
end
