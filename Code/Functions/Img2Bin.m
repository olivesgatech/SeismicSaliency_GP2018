function [ Th_FUT ] = Img2Bin( SS_UT, SS_Num)

    if (SS_Num <= 123)
        Thres_red  = 0.017;            % Threshold reduction 
    elseif (SS_Num <= 132)
        Thres_red  = 0.055;            % Threshold reduction 
    elseif (SS_Num <= 160)
        Thres_red  = 0.04;             % Threshold reduction 
    else
        Thres_red  = 0;
    end

    Th_FUT = im2bw(SS_UT, graythresh(SS_UT)-Thres_red);       
end

