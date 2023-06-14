% Load data
data = struct();

% Define file names and corresponding variable names
files = {'AZN.csv', 'JNJ.csv', 'MRNA.csv', 'PFE.csv'};
names = {'AstraZeneca', 'Johnson & Johnson', 'Moderna', 'Pfizer'};

% Load data into struct
for i = 1:length(files)
    tableData = readtable(files{i});
    data.(names{i}) = table2array(tableData(:, 5));
end

% Plot autocorrelation for each time series
figure;
for i = 1:length(names)
    subplot(2, 2, i);
    autocorr(data.(names{i}));
    title(names{i});
end

% Plot all time series
figure('Position', [100 100 800 800]);
g = gramm();

% Add data to the gramm object
for i = 1:length(names)
    g(1, i) = gramm('x', 1:length(data.(names{i})), 'y', data.(names{i}));
    g(1, i).axe_property('YLim', [30 190]);
    g(1, i).set_color_options('map', 'd3_10');
end

g.set_names('x', 'Day');
g.set_layout_options('legend', false);
g.draw();

% Test stationarity using Augmented Dickey Fuller
models = {'AR', 'TS'};
pValues = zeros(length(names), length(models));

for i = 1:length(names)
    for j = 1:length(models)
        [~, pValue] = adftest(data.(names{i}), 'model', models{j});
        pValues(i, j) = pValue;
    end
end

disp(pValues);  % Display p-values

% Calculate log returns
returns = struct();

for i = 1:length(names)
    returns.(names{i}) = diff(log(data.(names{i})));
end

% Plot ACF, PACF, and QQ plot for log returns
figure;
for i = 1:length(names)
    subplot(2, 2, i);
    autocorr(returns.(names{i}));
    title(names{i});
end

figure;
for i = 1:length(names)
    subplot(2, 2, i);
    parcorr(returns.(names{i}));
    title(names{i});
end

figure;
for i = 1:length(names)
    subplot(2, 2, i);
    qqplot(returns.(names{i}));
    title(names{i});
end

% Plot first differences of log returns
figure('Position', [100 100 800 800]);
g = gramm();

% Add data to the gramm object
for i = 1:length(names)
    g(1, i) = gramm('x', 1:length(returns.(names{i})), 'y', returns.(names{i}));
    g(1, i).axe_property('YLim', [-0.2 0.2]);
    g(1, i).set_color_options('map', 'd3_10');
end

g.set_names('x', 'Day');
g.set_layout_options('legend', false);
g.draw();

% Plot histograms of log returns
figure('Position', [100 100 800 600]);
g2 = gramm();

% Add histograms to the gramm object
for i = 1:length(names)
    g2(i) = gramm('x', returns.(names{i}));
    g2(i).stat_bin('normalization', 'probability', 'geom', 'overlaid_bar');
    g2(i).set_title(names{i});
    g2(i).set_color_options('map', 'brewer3');
end

g2.draw();
% Test stationarity using Augmented Dickey Fuller
pValues = zeros(length(names), 2);

for i = 1:length(names)
    [~, pValues(i, 1)] = adftest(returns.(names{i}));
    [~, pValues(i, 2)] = adftest(returns.(names{i}), 'model', 'TS');
end

pValues

% Function to plot observed and predicted values
plotObservedAndPredicted(data, prediction, model;

% Function to forecast future values
forecastFutureValues(model, data);

% Function to determine ARMA order
determineARMAOrder(data);

% ARMA order for AZN
determineARMAOrder(AZN_pct);

mod = arima(1, 0, 1);
AZN_mod = estimate(mod, AZN_pct);
AZN_res = infer(AZN_mod, AZN_pct);
AZN_pred = AZN_pct - AZN_res;

figure;
plotObservedAndPredicted(AZN_pct, AZN_pred, AZN_mod);
forecastFutureValues(AZN_mod, AZN_pct);

% ARMA order for JNJ
determineARMAOrder(JNJ_pct);

mod = arima(3, 0, 3);
JNJ_mod = estimate(mod, JNJ_pct);
JNJ_res = infer(JNJ_mod, JNJ_pct);
JNJ_pred = JNJ_pct - JNJ_res;

figure;
plotObservedAndPredicted(JNJ_pct, JNJ_pred, JNJ_mod);
forecastFutureValues(JNJ_mod, JNJ_pct);

% ARMA order for MRNA
determineARMAOrder(MRNA_pct);

mod = arima(2, 0, 2);
MRNA_mod = estimate(mod, MRNA_pct);
MRNA_res = infer(MRNA_mod, MRNA_pct);
MRNA_pred = MRNA_pct - MRNA_res;

figure;
plotObservedAndPredicted(MRNA_pct, MRNA_pred, MRNA_mod);
forecastFutureValues(MRNA_mod, MRNA_pct);

% ARMA order for PFE
determineARMAOrder(PFE_pct);

mod = arima(1, 0, 1);
PFE_mod = estimate(mod, PFE_pct);
PFE_res = infer(PFE_mod, PFE_pct);
PFE_pred = PFE_pct - PFE_res;

figure;
plotObservedAndPredicted(PFE_pct, PFE_pred, PFE_mod);
forecastFutureValues(PFE_mod, PFE_pct);

% Function to perform mean model diagnostic
performMeanModelDiagnostic(residuals, model);


% Mean model diagnostic for each company
companies = {'AZN', 'JNJ', 'MRNA', 'PFE'};
residuals = {AZN_res, JNJ_res, MRNA_res, PFE_res};
modModels = {AZN_mod, JNJ_mod, MRNA_mod, PFE_mod};

for i = 1:numel(companies)
    meanModelDiagnostic(residuals{i}, modModels{i});
end

% Function to plot autocorrelation and partial autocorrelation of squared residuals
plotResiduals(residuals, company);

% Plot autocorrelation and partial autocorrelation of squared residuals for each company
for i = 1:numel(companies)
    plotResiduals(residuals{i}, companies{i});
end

% GARCH and EGARCH model estimation for each company
pValues = {[0:4], [1:4], [0:4], [0:9]};
qValues = {[1:4], [1:4], [1:4], [1:9]};
arimaOrders = {[1, 0, 0], [3, 0, 3], [2, 0, 2], [4, 0, 4]};

for i = 1:numel(companies)
    % Plot autocorrelation and partial autocorrelation of squared residuals
    plotResiduals(residuals{i});

    % GARCH model estimation
    Mdl = arima(arimaOrders{i});
    CVarMdl = garch(pValues{i}, qValues{i});
    Mdl.Variance = CVarMdl;
    [EstMdl, ~, logL] = estimate(Mdl, eval([companies{i}, '_pct']), 'Display', 'off');

    % AIC and BIC calculation
    [aic, bic] = aicbic(logL, cellfun(@(p, q) p + q + 1, pValues{i}, qValues{i}), numel(eval([companies{i}, '_pct'])));

    % Display AIC and BIC matrices
    disp([companies{i}, ' AIC Matrix:']);
    disp(reshape(aic, numel(pValues{i}), numel(qValues{i})));
    disp([companies{i}, ' BIC Matrix:']);
    disp(reshape(bic, numel(pValues{i}), numel(qValues{i})));

    % Estimate GARCH model
    dy = eval([companies{i}, '_pct']);
    Mdl.Variance = garch(0, max(qValues{i}));
    GARCH = estimate(Mdl, dy);
}

% EGARCH model estimation
egarchModels = {[1, 1], [3, 3], [2, 2], [4, 4]};
for i = 1:numel(companies)
    Mdl = arima([1, 0, 0]);
    CVarMdl = egarch(egarchModels{i}(1), egarchModels{i}(2));
    Mdl.Variance = CVarMdl;
    [~, ~, logL] = estimate(Mdl, eval([companies{i}, '_pct']));
    disp(['EGARCH_', companies{i}, ' = ', num2str(logL)]);
end

% Forecasting
forecastData = {'AZN_future.csv', 'JNJ_future.csv', 'MRNA_future.csv', 'PFE_future.csv'};
forecastVariables = {'AZNf', 'JNJf', 'MRNAf', 'PFEf'};
forecastVariablePcts = {'AZNf_pct', 'JNJf_pct', 'MRNAf_pct', 'PFEf_pct'};

for i = 1:numel(companies)
    % Load forecast data
    forecastTable = readtable(forecastData{i});
    forecastArray = table2array(forecastTable(:, 5));
    eval([forecastVariables{i}, ' = forecastArray;']);

    % Calculate percentage change
    forecastLog = log(forecastArray);
    laggedForecast = lagmatrix(forecastLog, 1);
    forecastDiff = forecastLog - laggedForecast;
    eval([forecastVariablePcts{i}, ' = forecastDiff(2:end);']);
end

% Forecasting GARCH vs EGARCH
forecastCompanies = {'AZN', 'JNJ', 'MRNA', 'PFE'};
forecastVariables = {'AZNf_pct', 'JNJf_pct', 'MRNAf_pct', 'PFEf_pct'};
forecastGARCH = {'AG', 'JG', 'MG', 'PG'};
forecastEGARCH = {'AE', 'JE', 'ME', 'PE'};
forecastMSE = {'MSE_AG', 'MSE_AE', 'MSE_JG', 'MSE_JE', 'MSE_MG', 'MSE_ME', 'MSE_PG', 'MSE_PE'};
forecastMAE = {'MAE_AG', 'MAE_AE', 'MAE_JG', 'MAE_JE', 'MAE_MG', 'MAE_ME', 'MAE_PG', 'MAE_PE'};

for i = 1:numel(forecastCompanies)
    garchModel = GARCH;
    egarchModel = eval(['EGARCH_', forecastCompanies{i}]);
    pctData = eval([forecastVariables{i}, '(1:10)']);

    % GARCH forecast
    [dyG, dyMSE] = forecast(garchModel, 10, pctData);
    s = dyG - pctData;
    eval([forecastMSE{2*i-1}, ' = sum(s.^2) / 10;']);
    eval([forecastMAE{2*i-1}, ' = sum(abs(s)) / 10;']);

    % EGARCH forecast
    [dyE, dyMSE] = forecast(egarchModel, 10, pctData);
    s = dyE - pctData;
    eval([forecastMSE{2*i}, ' = sum(s.^2) / 10;']);
    eval([forecastMAE{2*i}, ' = sum(abs(s)) / 10;']);
end

% Correlation
correlationCompanies = {'AZN_JNJ', 'AZN_MRNA', 'AZN_PFE', 'PFE_JNJ', 'MRNA_JNJ', 'PFE_MRNA'};
se = cell(1, 6);
CI = cell(1, 6);

figure;
for i = 1:numel(correlationCompanies)
    subplot(3, 2, i);
    correlationData = eval([correlationCompanies{i}, '_pct']);
    corrMatrix = corr(correlationData);
    [se{i}, CI{i}] = bootstrap_R2(correlationData);
    histogram(se{i});
    title(correlationCompanies{i});
    CI{i}
end

figure;
g6 = gramm;
for i = 1:numel(correlationCompanies)
    g6(i) = gramm('x', se{i});
    g6(i).stat_bin('fill', 'face');
    g6(i).set_title(correlationCompanies{i});
    g6(i).set_names('x', 'Standard Error', 'y', 'Frequency');
end
g6.set_title('Standard Error Distribution');
g6.set_color_options('map', 'd3_10');
g6.draw();
