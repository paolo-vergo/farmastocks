function performMeanModelDiagnostic(data)
    % Function to perform diagnostic analysis for the mean model

    % Fit the mean model (ARMA(0,0))
    model = fitARMA(data, [0, 0]);

    % Perform diagnostic analysis
    residuals = model.Residuals.Raw;

    % Display summary statistics of the residuals
    disp('Summary Statistics:');
    disp(summary(residuals));

    % Plot the residuals
    plotResiduals(residuals);
end
