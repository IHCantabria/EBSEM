function metVal = Objective(model, settings, params, MetObj, ENS)

if model == "Y09"
    Y = yates09(settings, exp(params));
elseif model == "MD"
    Y = millerDean04(settings, exp(params));
elseif model == "SF"
    Y = shorefor(settings, exp(params));
end

YYsl = Y(ENS.indexes);

if MetObj == "Pearson"
    metVal = 1 -  abs(sum((YYsl-mean(YYsl)).*(ENS.Yobs - mean(ENS.Yobs)))/...
        (sqrt(sum((YYsl-mean(YYsl)).^2).*sum((ENS.Yobs-mean(ENS.Yobs)).^2))));
elseif MetObj == "RMSE"
    metVal = sqrt(mean((YYsl - ENS.Yobs).^2));
elseif MetObj == "MSS"
    metVal = sum((YYsl - ENS.Yobs).^2)/length(YYsl)/(var(YYsl)+var(ENS.Yobs)+(mean(YYsl)-mean(ENS.Yobs))^2);
elseif MetObj == "BSS"
    metVal = (mean((YYsl - ENS.Yobs).^2) - mean((YYref - ENS.Yobs).^2))/mean((YYref - ENS.Yobs).^2);
end