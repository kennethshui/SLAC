function [mask, RfacR, RfacF, RFD]= OSS_revised_complex_gs(F2D,supp_y, supp_x, iter,beta,showim,modelimg,hiofirst) 
%function [mask, RfacR, RfacF, RFD]= OSS_revised_complex_gs (F2D,supp_y, supp_x, iter,beta,showim,modelimg,hiofirst)
%OSS Source code modified on 07/13/2022, based on 'OSS_revised_complex.m' on 09/09/22, to add real-space magnitude constraint (GS algorithm).
%F2D:           Diffraction pattern of the oversampled vesicle image, with the central speckle removed, and 5% Poisson noise (Rnoise) added as defined above.
%supp_y [52]:   Size of rectangular support along the first (y) dimension, always placed at the center of image.
%supp_x [28]:   Size of rectangular support along the second (x) dimension.
%iter [2000]:   Total number of iterations performed by the OSS reconstruction.
%beta [0.9]:    Parameter used in a typical HIO reconstruction (see Fienup, Opt. Let. 1978).
%showim [1]:    Defines whether or not images showing progress are displayed (1 is shown, 0 not shown).
%modelimg:      Image of the oversampled vesicle model used for the OSS reconstruction.
%hiofirst:      =0: use OSS. =1: use HIO method.
%Output:
%mask:          An image of the mask used as the support in the OSS reconstruction, it has the same size as F2D.
%RfacR:         Error measurement in real space, defined above, computed alongside Rf, 1000x1 vector.
%RfacF:         Error measurement in Fourier space, defined above, computed every 2 iterations, 1000x1 vector.
%RFD:           Final result image obtained from the OSS reconstruction.
%Usage:         load model.mat; load Ten_perc_noise_Diffraction_Pattern.mat;
%               [mask, RfacR, RfacF, RFD]= OSS_samplecode (diffpat,60,30,2000,0.9,1,model,0);

close all;

%% General assignment of variables
[Rsize,Csize] = size(F2D)
R2D=zeros(Rsize,Csize,10,'single');    %10 frame of 2-D real sapce data.
toperrs=single(1:10:100);
kfilter=zeros(Rsize,Csize,'single');
realgroup=zeros(Rsize,Csize,2,'single');
realgroup(:,:,1)=modelimg;
 
%% Assign variables
%stopper = find(F2D==-1);   %Magnitude of K-space should be non-negative. Negative number introduced could be due to additive noise.
filtercount=10*1;
filtnum=0;
store=0;
 
%% Define support
Rcenter = ceil(Rsize/2);
Ccenter = ceil(Csize/2);
Rsupport = supp_y;
Csupport = supp_x;
half_Rsupport = ceil(Rsupport/2);
half_Csupport = ceil(Csupport/2);
support = zeros(Rsize,Csize,'single');
%support(Rcenter-half_Rsupport+1:Rcenter+half_Rsupport-1,Ccenter-half_Csupport+1:Ccenter+half_Csupport-1) = 1;  %Original support code.
support(Rcenter-half_Rsupport+1:Rcenter+half_Rsupport,Ccenter-half_Csupport+1:Ccenter+half_Csupport) = 1;  %Modified support code.
mask=support;

%% Draw the real space image with support region.
figure(3)
imagesc((abs(modelimg)));
title('Magnitude of input complex data');
figure(4)
imagesc(angle(modelimg));
grid on;
title('Phase of input complex data');

 
%% Compute filter parameter alpha 
X=1:iter; 
FX=(filtercount+1-ceil(X*filtercount/iter))*ceil(iter/(1*filtercount)); 
FX=((FX-ceil(iter/filtercount))*(2*Rsize)/max(FX))+(2*Rsize/10);
figure(98), plot(X,FX), axis tight; title('OSS Filter Size v. Iterations'); grid on;

 
%% Generate initial filter
for kk=1:Rsize
    for jj=1:Csize
        kfilter(kk,jj)=exp( -( ( (sqrt((kk-Rcenter)^2+(jj-Ccenter)^2).^2) ) ./ (2* FX(1)^2) ));
    end
end
kfilter=kfilter/max(max(kfilter));
 
%% Assign random phases
%rng('shuffle','twister');
%rng(1, 'twister');
phase_angle = rand(Rsize,Csize,'single'); 
%phase_angle = -0.117*ones(Rsize, Csize);  %Don't use random phase!! Use 0.1 rad is a good starting point!
%phase_angle = zeros(Rsize, Csize);
%phase_angle = randn(Rsize, Csize, 'single');
 
%% Define initial k, r space
initial_k=F2D; %initial_k(initial_k==-1)=0;
k_space = initial_k.*exp(1i*phase_angle);
buffer_r_space = single((ifftn(ifftshift(k_space))));    %Remove the real number constraint.
 
%% Preallocate error arrays
RfacF = zeros(ceil(iter/2),1,'single');  counter1=0; errorF=1;
RfacR = zeros(ceil(iter/2),1,'single');  counter2=0; errorR=1;
 
%% Image display argument
if showim==1    
    figure(1),
end
if nargin<7
    hiofirst=0;
end
 
%% OSS iterations
for iteration = 1:iter
    
        %% OSS without Support & Positivity constraint
        r1_space = (ifftn(ifftshift(k_space)));
        
        %Impose real space magnitude constraint
        r_space = abs(modelimg).*exp(1j*angle(r1_space)); 
        sample = r_space.*mask;    %Support constraint.
        
        r_space = buffer_r_space-beta*r_space;
        %sample(sample<0)=r_space(sample<0);  %Positivity constraint.
        %sample(imag(sample)<=0) = r_space(imag(sample)<=0);  %Impose imaginary part of real-space be positive constraint.
        
        %% Apply frequency filter (OSS)
        if hiofirst==0 || iteration>ceil(iter/filtercount)
            for kk=1:Rsize
                for jj=1:Csize
                    kfilter(kk,jj)=exp( -( ( (sqrt((kk-Rcenter)^2+(jj-Ccenter)^2).^2) ) ./ (2* FX(iteration)^2) ));
                end
            end
            kfilter=kfilter/max(max(kfilter));
            ktemp=fftshift(fftn(r_space));
            ktemp=ktemp.*kfilter;
            r_space=single((ifftn(ifftshift(ktemp))));
        end
 
        %% Use best result from last filter
        if mod(iteration,ceil(iter/filtercount))==0   %At 200, 400, 600, ..., 2000, total 10 set of data
            r_space=R2D(:,:,filtnum);       %Note: initial vallue of filtnum=0. At 200, 400,... iterations, save filtnum R2D data to r_space.
        else
            r_space(mask==1)=sample(mask==1);  %If not at 200, 400, ..., copy sample data of support region to r_space support region. 
        end

        
        %% Update reconstruction
        buffer_r_space = r_space;
        k_space = fftshift(fftn(r_space));  
        phase_angle = angle(k_space);
  
        %stopper_k_space = k_space(stopper);   %??  
        k_space = F2D.*exp(1i*phase_angle);   %Impose K-space magnitude here.
        %k_space(stopper) = stopper_k_space;   %??
                   
        %% Calculate errors, computed every 2 iterations   
        if rem(iteration,2)==0  
            
            %% Calculate error in reciprocal(Fourier) space
            Ktemp = sample;
            Ktemp = abs(fftshift(fftn(Ktemp)));
            errorF = sum(sum(abs(Ktemp(F2D~=-1)-F2D(F2D~=-1)))) / sum(sum(abs(F2D(F2D~=-1))));    %F2D is the input diffraction pattern magnitude.
            counter1=counter1+1; RfacF(counter1) = errorF;
            
            %% Determine interations with best error
            filtnum=ceil(iteration*filtercount/iter);
            %toperrs(filtnum)
            %errorF
            %store
            if errorF<= toperrs(filtnum) && iteration>store+2
                toperrs(filtnum)=errorF;
                R2D(:,:,filtnum)=r_space;
                store=iteration;
            end
            %pause;
            %% Calculate error in real space    
            realgroup(:,:,2)=sample;  %Original code.
            %realgroup(:,:,2) = r_space;
            realgroup2=realgroup(Rcenter-half_Rsupport-1:Rcenter+half_Rsupport+1,Ccenter-half_Csupport-1:Ccenter+half_Csupport+1,:);
            realgroup2 = align2(realgroup2,0,iteration);
            %[realgroup2]=align2(realgroup2,1,iteration);
            %pause;
            errorR = sum(sum(abs(realgroup2(:,:,1)-realgroup2(:,:,2)))) / sum(sum(realgroup2(:,:,1)));
            counter2=counter2+1; RfacR(counter2) = errorR;
            
            %% Figure shows progress
            if showim==1
                figure (1),
                if hiofirst==0
                    subplot(2,2,1), imagesc(abs(squeeze(realgroup2(:,:,1)))), axis image, title(strcat(int2str(FX(iteration)),'--OSS'));
                else
                    subplot(2,2,1), imagesc(abs(squeeze(realgroup2(:,:,1)))), axis image, title('HIO');
                    %subplot(2,2,1), imagesc(abs(squeeze(modelimg))), axis image, title('HIO');
                end
                subplot(2,2,2), imagesc(abs(squeeze(realgroup2(:,:,2)))), axis image, title(strcat('Iteration = ', int2str(iteration)));
                subplot(2,2,3), plot(RfacF), axis([0 ceil(iteration/2) 0 0.8]), title(strcat('F error=',int2str(errorF*100),'%')); grid on;
                subplot(2,2,4), plot(RfacR), axis([0 ceil(iteration/2) 0 1.0]), title(strcat('R error=',int2str(errorR*100),'%')); grid on;
                drawnow
            end
        end
end
 
%% Save results
if rem(iteration,iter)==0
    save ('MyR2D.mat','R2D','RfacF','RfacR','mask','toperrs','r_space');
end
 
%% Show image: sum of best 4 steps
 s=find(toperrs==min(toperrs));
 RFD= squeeze(R2D(:,:,s)); 
%RFD = squeeze(R2D(:,:,filtercount));  %No need to sort. Always choose the last one.

if showim==1
     figure(2),
     imagesc(abs(RFD));
     title('Retrieved image magnitude')
    subplot(2,2,1), imagesc(abs(RFD)), axis image, title('Retrieved image magnitude at Real Space'); 
    subplot(2,2,2), imagesc(squeeze(mask)), axis image, title('Mask/support area');
    subplot(2,2,3), plot(RfacF), axis tight, title('F-space error'); grid on;
    %subplot(2,2,4), plot(RfacR), axis tight, title('R-space error'); grid on;
    subplot(2,2,4), imagesc(angle(RFD)), axis image, title('Retrived image phase at Real Space');
end
