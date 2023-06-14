function plotResiduals(residuals)
    % Function to plot the residuals of a time series model

    % Create a time axis for the residuals
    timeAxis = 1:numel(residuals);

    % Plot the residuals over time
    plot(timeAxis, residuals);

    % Set plot title and labels
    title('Residuals Plot');
    xlabel('Time');
    ylabel('Residuals');
end
