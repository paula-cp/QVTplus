function [mdl,plotData] = fit_linXYData(xData,yData)
        mdl = fitlm(xData,yData,'linear');
        rng=[min(xData):(max(xData)-min(xData))/50:max(xData)]';
        [Ymain,yci] = predict(mdl,rng);
        plotData=[rng,Ymain,yci];
end