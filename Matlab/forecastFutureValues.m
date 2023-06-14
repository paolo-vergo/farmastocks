function forecastFutureValues(data, order)
    % Function to forecast future values using ARMA model

    % Perform ARMA model fitting
    model = fitARMA(data, order);

    % Forecast future values
    futureValues = forecast(model, numel(data));

    % Display the forecasted values
    disp(futureValues);
end
