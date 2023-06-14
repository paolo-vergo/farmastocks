function order = determineARMAOrder(data)
    % Function to determine the order of the ARMA model

    % Perform model order selection using AIC criterion
    orders = [1, 1]; % Initialize the minimum order as ARMA(1,1)
    minAIC = inf;

    for p = 1:3
        for q = 1:3
            try
                model = fitARMA(data, [p, q]);
                AIC = model.ModelCriterion.AIC;
                if AIC < minAIC
                    minAIC = AIC;
                    orders = [p, q];
                end
            catch
                continue; % Skip invalid orders
            end
        end
    end

    order = orders;
end
