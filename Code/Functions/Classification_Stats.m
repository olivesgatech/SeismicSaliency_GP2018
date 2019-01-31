function [ Accuracy_Sal,Precision_Sal,Recall_Sal,F_score,G_score] = Classification_Stats( Salt_Dome_sal, Ground_Truth )

    [ TP, FP, TN, FN ]   = ROC_calc( Salt_Dome_sal, Ground_Truth);

    True_Pos_Rate     = TP / (TP + FN);   % Sensitivity, Recall
    False_Neg_Rate    = FN / (TP + FN);   % Miss rate, 1-True_Pos_Rate
    False_Pos_Rate    = FP / (TN + FP);   % Fall out,  1-True_neg_Rate
    True_Neg_Rate     = TN / (TN + FP);   % Specificity,

    Accuracy_Sal      = (TP + TN) / (TP + FP + FN + TN);
    Precision_Sal     = TP / (TP + FP);   % Positive Predictive Value
    Recall_Sal        = True_Pos_Rate;   % Recall
    
    Neg_Predictive_val   = TN / (TN + FN);
    False_Discovery_Rate = FP / (TP + FP);    % 1-Positive Predictive Value
    False_Omission_Rate  = TP / (TP + FP);

    Pos_Likelihood_Ratio = True_Pos_Rate / False_Pos_Rate;
    Neg_Likelihood_Ratio = False_Neg_Rate / True_Neg_Rate;

    Diagnostic_Odd_Ratio = Pos_Likelihood_Ratio / Neg_Likelihood_Ratio;

    F1_score = (2*TP) / ((2*TP) + FP + FN);
    F_score  = (2*Precision_Sal*True_Pos_Rate) / (Precision_Sal + True_Pos_Rate);

    G_score  = sqrt(Precision_Sal*Recall_Sal);
    
end

