function sim_curve = SalSIM_cmpr_mod_v01(Curve_1, Curve_2 )

%% Function start
[rows, cols] = size(Curve_1);

Curve_1(end,:) = 0;
Curve_2(end,:) = 0;

%% Get Points
% Trace all points of the salt-dome boundary in "gTruth"
pts_Curve_1 = pointOrder(Curve_1);
% Trace all points of the salt-dome boundary in "result"
pts_Curve_2 = pointOrder(Curve_2);

if (length(pts_Curve_1) < length(pts_Curve_2))
    % Swap curves
    Temp_crv = Curve_1;
    Curve_1  = Curve_2;
    Curve_2  = Temp_crv;
    pts_Curve_1 = pointOrder(Curve_1);
    pts_Curve_2 = pointOrder(Curve_2);
end

%% Distance Calculation

Segments_dist = [];
% Compare the segments of salt-dome boundaries
len = 20; % Number of points in each segment of pts_gt

for i = 1+round(len/2):len/2:length(pts_Curve_1)-round(len/2)

    seg1 = pts_Curve_1(i-round(len/2):i+round(len/2),:);

    % According to the beginning and ending points of seg1 to determine the
    % beginning and ending points of seg2
    [~, seg2_beginIdx] = min(sum(abs(pts_Curve_2 - repmat(seg1(1,:),   length(pts_Curve_2), 1)), 2));
    [~, seg2_endIdx]   = min(sum(abs(pts_Curve_2 - repmat(seg1(end,:), length(pts_Curve_2), 1)), 2));

    maxIdx = max(seg2_beginIdx, seg2_endIdx);
    minIdx = min(seg2_beginIdx, seg2_endIdx);
    seg2 = pts_Curve_2(minIdx:maxIdx,:);
    
    % Calculate the Frechet Distance between seg1 and seg2
    dist_Fr = DiscreteFrechetDist(seg1, seg2);
    Segments_dist = [Segments_dist; dist_Fr];
end

%%

% Distance between the beginning points of "gTruth" and "Detected Bound"
[Row_Curve_1, Col_Curve_1] = ind2sub([rows, cols], find(Curve_1(:) == 1));
beginCol_Curve_1      = Col_Curve_1(1);
endCol_Curve_1        = Col_Curve_1(end);
beginCols_Curve_1     = find(Col_Curve_1 == beginCol_Curve_1);
beginRow_Curve_1      = Row_Curve_1(beginCols_Curve_1(end));
endRow_Curve_1        = Row_Curve_1(end);
lengthCol_Curve_1     = endCol_Curve_1 - beginCol_Curve_1 + 1;

% Distance between the ending points of "gTruth" and "Detected Bound"
[Row_Curve_2, Col_Curve_2] = ind2sub([rows, cols], find(Curve_2(:) == 1));
beginCol_Curve_2      = Col_Curve_2(1);
endCol_Curve_2        = Col_Curve_2(end);
beginCols_Curve_2     = find(Col_Curve_2 == beginCol_Curve_2);
beginRow_Curve_2      = Row_Curve_2(beginCols_Curve_2(end));
endRow_Curve_2        = Row_Curve_2(end);
lengthCol_Curve_2     = endCol_Curve_2 - beginCol_Curve_2 + 1;

distBegin   = sqrt((beginRow_Curve_1 - beginRow_Curve_2).^2 + (beginCol_Curve_1 - beginCol_Curve_2).^2);
distEnd     = sqrt((endRow_Curve_1   -   endRow_Curve_2).^2 + (endCol_Curve_1   -   endCol_Curve_2).^2);

lambda_Begin  = abs((beginCol_Curve_1 - beginCol_Curve_2 + 1))/lengthCol_Curve_1;
lambda_End    = abs((endCol_Curve_2   - endCol_Curve_1   + 1))/lengthCol_Curve_1;

% Statistical results of the distances between segments
distMean  = mean(Segments_dist);
distStd   = std(Segments_dist);
distMax   = max(Segments_dist);
penalty   = distBegin*lambda_Begin + distEnd*lambda_End;

alpha = 0.006;
beta  = 0.0015;
idx1  = exp(-alpha*(distMean+distStd));
idx2  = exp(-beta*(distMax+penalty));
idx   = exp(-alpha*(distMean + distStd)-beta*(distMax + penalty));
distAll = distMax;

sim_curve = idx;





