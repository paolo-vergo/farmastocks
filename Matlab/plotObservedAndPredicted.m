function plotObservedAndPredicted(data, prediction, model)
    plot(data, '-b');
    hold on;
    plot(prediction, '-r');
    title('1-Step Ahead', 'interpreter', 'latex');
    legend('Observed', 'Forecast', 'Location', 'Best', 'interpreter', 'latex');
    grid on;
    set(gca, 'FontSize', 20);
end
