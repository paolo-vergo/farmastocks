function [se, CI] = bootstrap_R2(v1, v2)
    %INPUT
    % v1, v2 = vectors of equal length
    %
    %OUTPUT
    % se = vector of correlations obtained via bootstrap
    % CI = 95% confidence interval for correlation

    % Fit linear regression model
    mdl = fitlm(v2, v1);
    fitted = mdl.Fitted;
    residuals = table2array(mdl.Residuals(:, 1));
    R2 = mdl.Rsquared.Ordinary;

    % Bootstrap procedure
    numBootstraps = 10000;
    se = bootstrp(numBootstraps, @(bootr) corr(fitted + bootr, v2), residuals);

    % Calculate mean and standard deviation of bootstrap correlations
    mean_se = mean(se);
    std_se = std(se);

    % Calculate confidence interval
    alpha = 0.05;
    z = norminv(1 - alpha / 2);
    CI = [mean_se - z * std_se, mean_se, mean_se + z * std_se];
end
