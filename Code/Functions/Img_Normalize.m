function [ Norm_SS ] = Img_Normalize( grayImage , rangenorm )


originalMinValue = double(min(min(grayImage)));
originalMaxValue = double(max(max(grayImage)));
originalRange = originalMaxValue - originalMinValue;

% Get a double image in the range 0 to +1
desiredMin = rangenorm(1);
desiredMax = rangenorm(2);
desiredRange = desiredMax - desiredMin;
Norm_SS = desiredRange * (double(grayImage) - originalMinValue) / originalRange + desiredMin;

end

