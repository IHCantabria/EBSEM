function [RMSE, RP, MSS] = MultiObjective(model, settings, params, ENS)

if model == "Y09"
    Y = yates09(settings, exp(params));
elseif model == "MD"
    Y = millerDean04(settings, exp(params));
elseif model == "SF"
    Y = shorefor(settings, exp(params));
end

YYsl = Y(ENS.indexes);

RP = 1 -  abs(sum((YYsl-mean(YYsl)).*(ENS.Yobs - mean(ENS.Yobs)))/...
        (sqrt(sum((YYsl-mean(YYsl)).^2).*sum((ENS.Yobs-mean(ENS.Yobs)).^2))));
RMSE = sqrt(mean((YYsl - ENS.Yobs).^2))/5;
MSS = sum((YYsl - ENS.Yobs).^2)/length(YYsl)/(var(YYsl)+var(ENS.Yobs)+(mean(YYsl)-mean(ENS.Yobs))^2);