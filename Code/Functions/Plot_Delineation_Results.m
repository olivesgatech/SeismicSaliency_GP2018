function Plot_Delineation_Results( SS_Display_result, SS_dataset, crop_val, xi, yi, SD_2DGoT, SD_3DGoT, SD_Oslu, SD_Aqrawif, SD_Codebook, SD_Sal )

    figure
    clf
    Temp_disp = imcrop(SS_dataset(:,:,SS_Display_result),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
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
    hold on
    % 2D GoT
    Salt_Dome_2Dd    = imcrop(SD_2DGoT(:,:,SS_Display_result),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    contour(Salt_Dome_2Dd, [0.5 0.5],'c','LineWidth',2);
    % 3D GoT
    Salt_Dome_3Dd    = imcrop(SD_3DGoT(:,:,SS_Display_result),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    contour(Salt_Dome_3Dd, [0.5 0.5],'b','LineWidth',2);
    % Oslu
    Oslu_SDd = imcrop(SD_Oslu(:,:,SS_Display_result),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    contour(Oslu_SDd, [0.5 0.5],'y','LineWidth',2);
    % Aqrawi
    Aqrawi_SDd = imcrop(SD_Aqrawif(:,:,SS_Display_result),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    contour(Aqrawi_SDd, [0.5 0.5],'m','LineWidth',2);
    % Codebook
    Codebook_SDd =  imcrop(SD_Codebook(:,:,SS_Display_result),[ crop_val crop_val yi-2*crop_val xi-2*crop_val]);
    contour(Codebook_SDd, [0.5 0.5],'k','LineWidth',2);
    % Saliency
    contour(SD_Sal(:,:,SS_Display_result), [0.5 0.5],'r','LineWidth',2);
    % Ground Truth
    load(['Ground_Truth_SS', num2str(SS_Display_result)],'-mat')
    Ground_Truthd  = imcrop(GT_SS,[ crop_val crop_val yi-2*crop_val xi-2*crop_val]); 
    Ground_TruthBd = bwperim(imcrop(Ground_Truthd,[ crop_val crop_val yi-2*crop_val xi-2*crop_val]));
    Ground_TruthBd(end,:) = false;
    contour(Ground_Truthd,[0.5,0.5],'g','linewidth',2)
    hold off

end

