function [ measure ] = fitMeasure( YPred, YTest, YTrain )

% Tjur R2
measure = sum(YPred.*YTest)./sum(YTest) -  sum(YPred.*(~YTest))./sum(~YTest);

% Deviance
% measure  = mean(-2*log(YTest.*YPred + (1-YTest).*(1-YPred)));

end

