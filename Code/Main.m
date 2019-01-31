%% Saliency based SD Detection
% ===============================================================================
%  Please, cite the following papers in your publication if you use this script.
% ===============================================================================
% M. Shafiq, T. Alshawi, Z. Long, and G. AlRegib, “The role of visual saliency in the automation of seismic 
% interpretation,” Geophysical Prospecting, 2017.

% M. Shafiq, T. Alshawi, Z. Long, and G. AlRegib, "SalSi: A New Seismic Attribute For Salt Dome 
% Detection," IEEE Intl. Conf. on Acoustics, Speech and Signal Processing (ICASSP), Shanghai, China, Mar. 20-25, 2016.

% If you have found any bugs or have any questions/suggestions, 
% please contact amirshafiq@gatech.edu 

%% Init
clc
clear
close all

%% Add Path
addpath(genpath('Functions'))
addpath(genpath('Mat Files'))

%% Read Saliency Video
load salMapAllSlices

%% Read Seismic database
load  salt2_inline.mat;
SS_dataset = salt2_inline;
clear salt2_inline

load SD_Aqrawi_fixed
load SD_Oslu
load SD_2DGoT
load SD_3DGoT
load SD_Codebook

load Frechet_Sim_Updated
load CurveD_Sal_GP

%% Saliency based Salt Dome detection
Display_Overlay_Sal = 0;
Intermediate_res    = 1;
[xi, yi]            = size(SS_dataset(:,:,100));

%% Sim Parameters
SS_UT_num  = 104:160; 

%% Threshold Loop
for SS_loop = 1:length(SS_UT_num)   
    SS_Num = SS_UT_num(SS_loop);
    SS_UT  = salMapNormalized(:,:,SS_Num);
    IUT    = SS_dataset(:,:,SS_Num);
    
    crop_val   = 10;               % Crop Value
    Dil_b4_RG  = 3;                % Pre-processing
    Dil_Aft_RG = 12;               % Post-processing
    SS_UT   = imcrop(SS_UT,[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    IUT     = imcrop(IUT,[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);

    Th_FUT(:,:,SS_loop) = Img2Bin( SS_UT, SS_Num); % Binary conversion
    Sal_OP(:,:,SS_loop) = imclose(Th_FUT(:,:,SS_loop), strel('disk',Dil_b4_RG));   % Close Saliency Maps
    Th_Bin_Vol(:,:,SS_loop) = Sal_OP(:,:,SS_loop); % Bin Volume
    
end

%% 3D Region Growing
load iSeed_pts_dataset
seed_pt = round(iSeed{SS_Num});
Seed_pt = [ seed_pt(2),seed_pt(1),1];

[Salt_Dome_3DRG, ~] = RegionGrowing_mod(Th_Bin_Vol,1,Seed_pt);  
Salt_Dome_Dil3D     = imdilate(Salt_Dome_3DRG, strel('disk',Dil_Aft_RG));  

%% Inlines Loop
for SS_loop = 1:length(SS_UT_num)   

    SS_Num = SS_UT_num(SS_loop);
    SS_UT  = salMapNormalized(:,:,SS_Num);
    IUT    = SS_dataset(:,:,SS_Num);

    SS_UT   = imcrop(SS_UT,[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    IUT     = imcrop(IUT,[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    
    % Load Ground Truth
    load(['Ground_Truth_SS', num2str(SS_Num)],'-mat')
    Ground_Truth  = imcrop(GT_SS,[ crop_val crop_val yi-2*crop_val xi-2*crop_val]); 
    Ground_TruthB = bwperim(imcrop(Ground_Truth,[ crop_val crop_val yi-2*crop_val xi-2*crop_val]));
    Ground_TruthB(end,:) = false;

    Salt_Dome_sal         = Salt_Dome_Dil3D(:,:,SS_loop);
    Salt_Dome_salB        = bwperim(imcrop(Salt_Dome_sal,[ crop_val crop_val yi-2*crop_val xi-2*crop_val]));
    Salt_Dome_salB(end,:) = false;
    SD_Sal(:,:,SS_Num)    = Salt_Dome_sal;

    %% Overlaying Salient boundary on SS
    % 2D GoT
    Salt_Dome_2D    = imcrop(SD_2DGoT(:,:,SS_Num),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    Salt_Dome_2DB        = bwperim(Salt_Dome_2D);
    Salt_Dome_2DB(end,:) = false;

    % 3D GoT
    Salt_Dome_3D    = imcrop(SD_3DGoT(:,:,SS_Num),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    Salt_Dome_3DB        = bwperim(Salt_Dome_3D);
    Salt_Dome_3DB(end,:) = false;

    % Oslu
    Oslu_SD = imcrop(SD_Oslu(:,:,SS_Num),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    Oslu_SDB        = bwperim(Oslu_SD);
    Oslu_SDB(end,:) = false;

    % Aqrawi
    Aqrawi_SD = imcrop(SD_Aqrawif(:,:,SS_Num),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    Aqrawi_SDB        = bwperim(Aqrawi_SD);
    Aqrawi_SDB(end,:) = false;

    % Codebook
    Codebook_SD =  imcrop(SD_Codebook(:,:,SS_Num),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    Codebook_SDB        = bwperim(Codebook_SD);
    Codebook_SDB(end,:) = false;

    Frech_sim_all(SS_loop, 2) = Frechet_2D_GoT(SS_Num);
    Frech_sim_all(SS_loop, 3) = Frechet_3D_GoT(SS_Num);
    Frech_sim_all(SS_loop, 4) = Frechet_Oslu(SS_Num);
    Frech_sim_all(SS_loop, 5) = Frechet_Aqrawi_Fixed(SS_Num);
    Frech_sim_all(SS_loop, 6) = Frechet_Codebook(SS_Num);
    Frech_sim_all(SS_loop, 7) = SalSIM_cmpr_mod_v01(Ground_TruthB,Salt_Dome_salB);

    fprintf('Frec_Sim of SS#%d using 2D GoT is    ==> %f \n',SS_Num, Frech_sim_all(SS_loop, 2));
    fprintf('Frec_Sim of SS#%d using 3D GoT is    ==> %f \n',SS_Num, Frech_sim_all(SS_loop, 3));
    fprintf('Frec_Sim of SS#%d using Oslu is      ==> %f \n',SS_Num, Frech_sim_all(SS_loop, 4));
    fprintf('Frec_Sim of SS#%d using Aqrawi is    ==> %f \n',SS_Num, Frech_sim_all(SS_loop, 5));
    fprintf('Frec_Sim of SS#%d using Codebook is  ==> %f \n',SS_Num, Frech_sim_all(SS_loop, 6));
    fprintf('Frec_Sim of SS#%d using Saliency is  ==> %f \n\n',SS_Num, Frech_sim_all(SS_loop, 7));

    % CurveD
    % fprintf('Computing CurveD \n');
    % Saliency_SimIdx_Curvature(SS_loop,:) = Curvature_Similarity(Ground_Truth, Salt_Dome_sal);
    % Saliency_SimIdx_Shape(SS_loop,:)     = Shape_Similarity(Ground_Truth, Salt_Dome_sal);
    % CurveD_Saliency(SS_loop,:)           = 0.5*Saliency_SimIdx_Curvature(SS_loop,:) + 0.5*Saliency_SimIdx_Shape(SS_loop,:);
    
    % Classification Accuracy Calculations
    [Acc_Sal(SS_loop), Prec_Sal(SS_loop), Rec_Sal(SS_loop), F_score_Sal(SS_loop), ...
                            G_score_Sal(SS_loop) ] = Classification_Stats( Salt_Dome_sal, Ground_Truth );
    [Acc_2DGoT(SS_loop), Prec_2DGoT(SS_loop), Rec_2DGoT(SS_loop), F_score_2DGoT(SS_loop),...
                            G_score_2DGoT(SS_loop)] = Classification_Stats( Salt_Dome_2D, Ground_Truth );
    [Acc_3DGoT(SS_loop), Prec_3DGoT(SS_loop), Rec_3DGoT(SS_loop), F_score_3DGoT(SS_loop), ...
                            G_score_3DGoT(SS_loop)] = Classification_Stats( Salt_Dome_3D, Ground_Truth );
    [Acc_Aqrawi(SS_loop), Prec_Aqrawi(SS_loop), Rec_Aqrawi(SS_loop),  F_score_Aqrawi(SS_loop), ...
                            G_score_Aqrawi(SS_loop)] = Classification_Stats( Aqrawi_SD, Ground_Truth );
    [Acc_Oslu(SS_loop), Prec_Oslu(SS_loop), Rec_Oslu(SS_loop), F_score_Oslu(SS_loop), ...
                            G_score_Oslu(SS_loop)] = Classification_Stats( Oslu_SD, Ground_Truth );
    [Acc_Codebook(SS_loop), Prec_Codebook(SS_loop), Rec_Codebook(SS_loop), F_score_Codebook(SS_loop), ...
                            G_score_Codebook(SS_loop)] = Classification_Stats( Codebook_SD, Ground_Truth );
    
    pause(0.01)
end

    
%% Display delineation results
Plot_Delineation_Results( 111, SS_dataset, crop_val, xi, yi, SD_2DGoT, SD_3DGoT, SD_Oslu, SD_Aqrawif, SD_Codebook, SD_Sal );
Plot_Delineation_Results( 123, SS_dataset, crop_val, xi, yi, SD_2DGoT, SD_3DGoT, SD_Oslu, SD_Aqrawif, SD_Codebook, SD_Sal );
Plot_Delineation_Results( 141, SS_dataset, crop_val, xi, yi, SD_2DGoT, SD_3DGoT, SD_Oslu, SD_Aqrawif, SD_Codebook, SD_Sal );
Plot_Delineation_Results( 159, SS_dataset, crop_val, xi, yi, SD_2DGoT, SD_3DGoT, SD_Oslu, SD_Aqrawif, SD_Codebook, SD_Sal );


%% Plot Frech of SS
SS_disp_axis = SS_UT_num + 249;

figure
plot(SS_disp_axis, Frech_sim_all(:,2),'-oc','linewidth', 2)
hold on
plot(SS_disp_axis, Frech_sim_all(:,3),'-+b','linewidth', 2)
plot(SS_disp_axis, Frech_sim_all(:,4),'-py','linewidth', 2)
plot(SS_disp_axis, Frech_sim_all(:,5),'-xm','linewidth', 2)
plot(SS_disp_axis, Frech_sim_all(:,6),'-sk','linewidth', 2)
plot(SS_disp_axis, Frech_sim_all(:,7),'-*r','linewidth', 2)
hold off
legend('Wang et al. (2015)','Shafiq et al. (2015)','Berthelot et al. (2013)','Aqrawi et al. (2011)','Amin and Deriche (2016)','Proposed method')
axis tight
xlabel('Seismic inline sections', 'FontSize', 30 , 'fontname','Times New Roman');
ylabel('SalSIM', 'FontSize', 30 , 'fontname','Times New Roman');
set(gca,'fontsize',30,'fontname','Times New Roman')  % set axis font size
axis tight
ylim([0.6 1])
grid on
grid minor

%% Plot CurveD
MkrSize = 7;
LineWd  = 2;
figure
clf
hold on
plot(SS_disp_axis, GoT2D_SimIdx(SS_UT_num),   '-oc','linewidth',LineWd)
plot(SS_disp_axis, GoT3D_SimIdx(SS_UT_num),   '-+b','linewidth',LineWd)
plot(SS_disp_axis, Oslu_SimIdx(SS_UT_num),    '-py','linewidth',LineWd)
plot(SS_disp_axis, Aqrawi_SimIdx(SS_UT_num),  '-xm','linewidth',LineWd)
plot(SS_disp_axis, Codebook_SimIdx(SS_UT_num),'-sk','linewidth',LineWd)
plot(SS_disp_axis, SalMulti_SimIdx(SS_UT_num),'-*r','linewidth',LineWd)
hold off
legend('Wang et al. (2015)','Shafiq et al. (2015)','Berthelot et al. (2013)','Aqrawi et al. (2011)','Amin and Deriche (2016)','Proposed method')
xlabel('Seismic inline sections', 'FontSize', 30 , 'fontname','Times New Roman');
ylabel('CurveD', 'FontSize', 30 , 'fontname','Times New Roman');
set(gca,'fontsize',30,'fontname','Times New Roman')  % set axis font size
axis tight
ylim([0 85])
grid on
box on
grid minor

%% Plot Accu, prec, Recall

figure
plot(SS_disp_axis, Acc_2DGoT, '-oc', 'LineWidth',2);
hold on
plot(SS_disp_axis, Acc_3DGoT, '-+b', 'LineWidth',2);
plot(SS_disp_axis, Acc_Oslu, '-py', 'LineWidth',2);
plot(SS_disp_axis, Acc_Aqrawi, '-xm', 'LineWidth',2);
plot(SS_disp_axis, Acc_Codebook, '-sk', 'LineWidth',2);
plot(SS_disp_axis, Acc_Sal, '-*r', 'LineWidth',2);
hold off
legend('Wang et al. (2015)','Shafiq et al. (2015)','Berthelot et al. (2013)','Aqrawi et al. (2011)','Amin and Deriche (2016)','Proposed method')
xlabel('Seismic inline sections', 'FontSize', 30 , 'fontname','Times New Roman');
ylabel('Accuracy', 'FontSize', 30 , 'fontname','Times New Roman');
set(gca,'fontsize',30,'fontname','Times New Roman')  % set axis font size
axis tight
ylim([0.85 1])
grid on
grid minor

figure
plot(SS_disp_axis, Prec_2DGoT, '-oc', 'LineWidth',2);
hold on
plot(SS_disp_axis, Prec_3DGoT, '-+b', 'LineWidth',2);
plot(SS_disp_axis, Prec_Oslu, '-py', 'LineWidth',2);
plot(SS_disp_axis, Prec_Aqrawi, '-xm', 'LineWidth',2);
plot(SS_disp_axis, Prec_Codebook, '-sk', 'LineWidth',2);
plot(SS_disp_axis, Prec_Sal, '-*r', 'LineWidth',2);
hold off
legend('Wang et al. (2015)','Shafiq et al. (2015)','Berthelot et al. (2013)','Aqrawi et al. (2011)','Amin and Deriche (2016)','Proposed method')
xlabel('Seismic inline sections', 'FontSize', 30 , 'fontname','Times New Roman');
ylabel('Precision', 'FontSize', 30 , 'fontname','Times New Roman');
set(gca,'fontsize',30,'fontname','Times New Roman')  % set axis font size
axis tight
ylim([0.75 1])
grid on
grid minor

figure
plot(SS_disp_axis, F_score_2DGoT, '-oc', 'LineWidth',2);
hold on
plot(SS_disp_axis, F_score_3DGoT, '-+b', 'LineWidth',2);
plot(SS_disp_axis, F_score_Oslu, '-py', 'LineWidth',2);
plot(SS_disp_axis, F_score_Aqrawi, '-xm', 'LineWidth',2);
plot(SS_disp_axis, F_score_Codebook, '-sk', 'LineWidth',2);
plot(SS_disp_axis, F_score_Sal, '-*r', 'LineWidth',2);
hold off
legend('Wang et al. (2015)','Shafiq et al. (2015)','Berthelot et al. (2013)','Aqrawi et al. (2011)','Amin and Deriche (2016)','Proposed method')
xlabel('Seismic inline sections', 'FontSize', 30 , 'fontname','Times New Roman');
ylabel('F-score', 'FontSize', 30 , 'fontname','Times New Roman');
set(gca,'fontsize',30,'fontname','Times New Roman')  % set axis font size
axis tight
ylim([0.85 1])
grid on
grid minor

fprintf('\n\n                  2DGoT    3DGoT    Oslu   Aqrawi  Codebook  Saliency \n');
fprintf('Mean Accuracy:    %0.4f  %0.4f   %0.4f  %0.4f   %0.4f    %0.4f \n',mean(Acc_2DGoT), mean(Acc_3DGoT),  mean(Acc_Oslu), mean(Acc_Aqrawi), mean(Acc_Codebook), mean(Acc_Sal));
fprintf('Mean Precision:   %0.4f  %0.4f   %0.4f  %0.4f   %0.4f    %0.4f \n',mean(Prec_2DGoT), mean(Prec_3DGoT),  mean(Prec_Oslu), mean(Prec_Aqrawi), mean(Prec_Codebook), mean(Prec_Sal));
fprintf('Mean Recall:      %0.4f  %0.4f   %0.4f  %0.4f   %0.4f    %0.4f \n',mean(Rec_2DGoT), mean(Rec_3DGoT),  mean(Rec_Oslu), mean(Rec_Aqrawi), mean(Rec_Codebook), mean(Rec_Sal));
fprintf('Mean F-measure:   %0.4f  %0.4f   %0.4f  %0.4f   %0.4f    %0.4f \n',mean(F_score_2DGoT), mean(F_score_3DGoT),  mean(F_score_Oslu), mean(F_score_Aqrawi), mean(F_score_Codebook), mean(F_score_Sal));
fprintf('Mean G-measure:   %0.4f  %0.4f   %0.4f  %0.4f   %0.4f    %0.4f \n',mean(G_score_2DGoT), mean(G_score_3DGoT),  mean(G_score_Oslu), mean(G_score_Aqrawi), mean(G_score_Codebook), mean(G_score_Sal));

fprintf('\n                  2DGoT    3DGoT    Oslu   Aqrawi  Codebook  Saliency \n');
fprintf('Std Accuracy:     %0.4f  %0.4f   %0.4f  %0.4f   %0.4f    %0.4f \n',std(Acc_2DGoT), std(Acc_3DGoT),  std(Acc_Oslu), std(Acc_Aqrawi), std(Acc_Codebook), std(Acc_Sal));
fprintf('Std Precision:    %0.4f  %0.4f   %0.4f  %0.4f   %0.4f    %0.4f \n',std(Prec_2DGoT), std(Prec_3DGoT),  std(Prec_Oslu), std(Prec_Aqrawi), std(Prec_Codebook), std(Prec_Sal));
fprintf('Std Recall:       %0.4f  %0.4f   %0.4f  %0.4f   %0.4f    %0.4f \n',std(Rec_2DGoT), std(Rec_3DGoT),  std(Rec_Oslu), std(Rec_Aqrawi), std(Rec_Codebook), std(Rec_Sal));
fprintf('Std F-measure:    %0.4f  %0.4f   %0.4f  %0.4f   %0.4f    %0.4f \n',std(F_score_2DGoT), std(F_score_3DGoT),  std(F_score_Oslu), std(F_score_Aqrawi), std(F_score_Codebook), std(F_score_Sal));
fprintf('Std G-measure:    %0.4f  %0.4f   %0.4f  %0.4f   %0.4f    %0.4f \n',std(G_score_2DGoT), std(G_score_3DGoT),  std(G_score_Oslu), std(G_score_Aqrawi), std(G_score_Codebook), std(G_score_Sal));

%% SalSIM stats
SS_start = SS_UT_num(1);
SS_stop  = SS_UT_num(end);

Mean_2DGoT_SalSIM        = mean(Frechet_2D_GoT(SS_start:SS_stop));
Mean_3DGoT_SalSIM        = mean(Frechet_3D_GoT(SS_start:SS_stop));
Mean_Aqrawi_SalSIM 	     = mean(Frechet_Aqrawi_Fixed(SS_start:SS_stop));
Mean_Oslu_SalSIM         = mean(Frechet_Oslu(SS_start:SS_stop));
Mean_Codebook_SalSIM     = mean(Frechet_Codebook(SS_start:SS_stop));
Mean_Sal_SalSIM          = mean(Frech_sim_all(:, 7));

% Std Dev
Std_2DGoT_SalSIM         = std(Frechet_2D_GoT(SS_start:SS_stop));
Std_3DGoT_SalSIM         = std(Frechet_3D_GoT(SS_start:SS_stop));
Std_Aqrawi_SalSIM        = std(Frechet_Aqrawi_Fixed(SS_start:SS_stop));
Std_Oslu_SalSIM          = std(Frechet_Oslu(SS_start:SS_stop));
Std_Codebook_SalSIM      = std(Frechet_Codebook(SS_start:SS_stop));
Std_Sal_SalSIM           = std(Frech_sim_all(:, 7));

fprintf('\n****************    SalSIM Similarity Mean    ****************\n');
fprintf('2D GoT      Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_2DGoT_SalSIM );
fprintf('3D GoT      Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_3DGoT_SalSIM );
fprintf('Oslu        Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_Oslu_SalSIM );
fprintf('Aqrawi      Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_Aqrawi_SalSIM );
fprintf('Codebook    Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_Codebook_SalSIM );
fprintf('Saliency    Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_Sal_SalSIM );

fprintf('\n\n****************    SalSIM Similarity Std    ****************\n');
fprintf('2D GoT      Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_2DGoT_SalSIM );
fprintf('3D GoT      Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_3DGoT_SalSIM );
fprintf('Oslu        Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_Oslu_SalSIM );
fprintf('Aqrawi      Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_Aqrawi_SalSIM );
fprintf('Codebook    Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_Codebook_SalSIM );
fprintf('Saliency    Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_Sal_SalSIM );

%% CurveD stats
% Mean
Mean_2DGoT        = mean(GoT2D_SimIdx(SS_start:SS_stop));
Mean_3DGoT        = mean(GoT3D_SimIdx(SS_start:SS_stop));
Mean_Aqrawi 	  = mean(Aqrawi_SimIdx(SS_start:SS_stop));
Mean_Oslu         = mean(Oslu_SimIdx(SS_start:SS_stop));
Mean_Codebook     = mean(Codebook_SimIdx(SS_start:SS_stop));
Mean_Saliency     = mean(SalMulti_SimIdx(SS_start:SS_stop));

Std_2DGoT         = std(GoT2D_SimIdx(SS_start:SS_stop));
Std_3DGoT         = std(GoT3D_SimIdx(SS_start:SS_stop));
Std_Aqrawi        = std(Aqrawi_SimIdx(SS_start:SS_stop));
Std_Oslu          = std(Oslu_SimIdx(SS_start:SS_stop));
Std_Codebook      = std(Codebook_SimIdx(SS_start:SS_stop));
Std_Saliency      = std(SalMulti_SimIdx(SS_start:SS_stop));

fprintf('\n****************    Similarity Mean    ****************\n');
fprintf('2D GoT      Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_2DGoT );
fprintf('3D GoT      Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_3DGoT );
fprintf('Oslu        Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_Oslu );
fprintf('Aqrawi      Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_Aqrawi );
fprintf('Codebook    Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_Codebook );
fprintf('Saliency    Mean for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Mean_Saliency );

fprintf('\n\n****************    Similarity Std    ****************\n');
fprintf('2D GoT      Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_2DGoT );
fprintf('3D GoT      Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_3DGoT );
fprintf('Oslu        Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_Oslu );
fprintf('Aqrawi      Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_Aqrawi );
fprintf('Codebook    Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_Codebook );
fprintf('Saliency    Standard Dev for SS#%d to SS#%d is ==> %0.4f \n',SS_start,SS_stop, Std_Saliency );

%% Intermediate Results

if Intermediate_res

    Display_SS = 160;

    % Orignal Section
    figure
    Temp_disp = imcrop(SS_dataset(:,:,Display_SS),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    Temp_DNorm = Img_Normalize(Temp_disp,[-1 1]);
    imshow(Temp_DNorm,[])
    colorbar
    xlabel('Crosslines')
    ylabel('Time (ms)')
    axis on 
    X_ticks = get(gca,'XTick');
    set(gca,'xTickLabel',X_ticks + 699);
    Y_ticks =  get(gca,'YTick');
    set(gca,'yTickLabel',(Y_ticks+325)*4);
    set(gca,'fontsize',24,'fontname', 'Times New Roman')

    % Sal Map
    figure
    Temp_disp  = imcrop(salMapNormalized(:,:,Display_SS),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    Temp_DNorm = Img_Normalize(Temp_disp,[0 1]);
    imshow(Temp_DNorm,[])
    colormap jet 
    colorbar
    caxis([0.09, 0.65])
    xlabel('Crosslines')
    ylabel('Time (ms)')
    axis on 
    X_ticks = get(gca,'XTick');
    set(gca,'xTickLabel',X_ticks + 699);
    Y_ticks =  get(gca,'YTick');
    set(gca,'yTickLabel',(Y_ticks+325)*4);
    set(gca,'fontsize',24,'fontname', 'Times New Roman')

    % Binary
    SS_Ind_Disp = Display_SS - SS_start + 1;
    figure
    Temp_disp  = Th_FUT(:,:,SS_Ind_Disp);
    imshow(Temp_disp,[])
    colorbar
    xlabel('Crosslines')
    ylabel('Time (ms)')
    axis on 
    X_ticks = get(gca,'XTick');
    set(gca,'xTickLabel',X_ticks + 699);
    Y_ticks =  get(gca,'YTick');
    set(gca,'yTickLabel',(Y_ticks+325)*4);
    set(gca,'fontsize',24,'fontname', 'Times New Roman')

    % RG
    seed_pt = round(iSeed{Display_SS});
    Seed_pt = [ seed_pt(2),seed_pt(1)];
    BI = Th_FUT(:,:,SS_Ind_Disp);
    [SD_SS, ~] = RegionGrowing_mod(BI,1,Seed_pt);  

    figure
    imshow(SD_SS,[])
    colorbar
    xlabel('Crosslines')
    ylabel('Time (ms)')
    axis on 
    X_ticks = get(gca,'XTick');
    set(gca,'xTickLabel',X_ticks + 699);
    Y_ticks =  get(gca,'YTick');
    set(gca,'yTickLabel',(Y_ticks+325)*4);
    set(gca,'fontsize',24,'fontname', 'Times New Roman')
    hold on
    plot(seed_pt(1)-10,seed_pt(2)-10,'r.','markersize',30)
    hold off

    % Post Proc
    figure
    imshow(Salt_Dome_Dil3D(:,:,SS_Ind_Disp),[])
    colorbar
    xlabel('Crosslines')
    ylabel('Time (ms)')
    axis on 
    X_ticks = get(gca,'XTick');
    set(gca,'xTickLabel',X_ticks + 699);
    Y_ticks =  get(gca,'YTick');
    set(gca,'yTickLabel',(Y_ticks+325)*4);
    set(gca,'fontsize',24,'fontname', 'Times New Roman')
    
    % Delineation
    figure
    Temp_disp = imcrop(SS_dataset(:,:,Display_SS),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    Temp_DNorm = Img_Normalize(Temp_disp,[-1 1]);
    imshow(Temp_DNorm,[])
    colorbar
    hold on
    contour(SD_Sal(:,:,Display_SS), [0.5 0.5],'g','LineWidth',2);
    hold off
    xlabel('Crosslines')
    ylabel('Time (ms)')
    axis on 
    X_ticks = get(gca,'XTick');
    set(gca,'xTickLabel',X_ticks + 699);
    Y_ticks =  get(gca,'YTick');
    set(gca,'yTickLabel',(Y_ticks+325)*4);
    set(gca,'fontsize',24,'fontname', 'Times New Roman')
    
end

