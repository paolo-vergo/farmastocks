%% Load all the needed data!

warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
AZN = readtable('AZN.csv');
AZN = AZN(:,5);
AZN = table2array(AZN);
JNJ = readtable('JNJ.csv');
JNJ = JNJ(:,5);
JNJ = table2array(JNJ);
MRNA = readtable('MRNA.csv');
MRNA = MRNA(:,5);
MRNA = table2array(MRNA);
PFE = readtable('PFE.csv');
PFE = PFE(:,5);
PFE = table2array(PFE);

%Plot autocorrelation
figure();
subplot(2,2,1);
autocorr(AZN);
title('AstraZeneca');
subplot(2,2,2);
autocorr(JNJ);
title('Johnson & Johnson');
subplot(2,2,3);
autocorr(MRNA);
title('Moderna');
subplot(2,2,4);
autocorr(PFE);
title('Pfizer');


% Plot all the time series

clear g;
g(1,1)=gramm('x',1:252,'y',AZN);
g(1,1).geom_line();
g(1,1).set_names('x','Day','y','AZN');
g(1,1).axe_property('YLim',[30 190]);
g(1,1).set_title('AstraZeneca');
g(1,2)=gramm('x',1:252,'y',JNJ);
g(1,2).geom_line();
g(1,2).axe_property('YLim',[30 190]);
g(1,2).set_names('x','Day','y','JNJ');
g(1,2).set_title('Johnson & Johnson');
g(2,1)=gramm('x',1:252,'y',MRNA);
g(2,1).geom_line();
g(2,1).axe_property('YLim',[30 190]);
g(2,1).set_names('x','Day','y','MRNA');
g(2,1).set_title('Moderna');
g(2,2)=gramm('x',1:252,'y',PFE);
g(2,2).geom_line();
g(2,2).set_names('x','Day','y','PFE');
g(2,2).set_title('Pfizer');
g(2,2).axe_property('YLim',[30 190]);
figure('Position',[100 100 800 800]);
g(1,1).set_color_options('map','d3_10');
g(1,2).set_color_options('map','brewer3');
g(2,1).set_color_options('map','winter');
g(2,2).set_color_options('map','brewer_pastel');
g.draw();

% TEST STATIONARITY
% Exploit the Augmented Dickey Fuller to check for the presence of unit root. 
% Overall we observe high p-values; therefore we cannot reject the null hypothesis 
% of unit root. Stationarity is not guaranteed.

%No TS
[~,pValue] = adftest(AZN); 
pValue
[~,pValue] = adftest(JNJ); 
pValue
[~,pValue] = adftest(MRNA); 
pValue
[~,pValue] = adftest(PFE); 
pValue

%TS
[~,pValue] = adftest(AZN,'model','TS'); 
pValue
[~,pValue] = adftest(JNJ,'model','TS'); 
pValue
[~,pValue] = adftest(MRNA,'model','TS'); 
pValue
[~,pValue] = adftest(PFE,'model','TS'); 
pValue

%% LOG RETURNS

AZN=log(AZN);
laggedAZN = lagmatrix(AZN,1);
diff=AZN-laggedAZN; 
AZN_pct=diff;
AZN_pct(1,:)=[];

JNJ=log(JNJ);
laggedJNJ = lagmatrix(JNJ,1);
diff=JNJ-laggedJNJ;
JNJ_pct=diff;
JNJ_pct(1,:)=[];

MRNA=log(MRNA);
laggedMRNA = lagmatrix(MRNA,1);
diff=MRNA-laggedMRNA;
MRNA_pct=diff;
MRNA_pct(1,:)=[];

PFE=log(PFE);
laggedPFE = lagmatrix(PFE,1);
diff=PFE-laggedPFE;
PFE_pct=diff;
PFE_pct(1,:)=[];

%plot ACF and PACF and QQ_plot

figure();
subplot(2,2,1);
autocorr(AZN_pct);
title('AstraZeneca');
subplot(2,2,2);
autocorr(JNJ_pct);
title('Johnson & Johnson');
subplot(2,2,3);
autocorr(MRNA_pct);
title('Moderna');
subplot(2,2,4);
autocorr(PFE_pct);
title('Pfizer');

figure();
subplot(2,2,1);
parcorr(AZN_pct);
title('AstraZeneca');
subplot(2,2,2);
parcorr(JNJ_pct);
title('Johnson & Johnson');
subplot(2,2,3);
parcorr(MRNA_pct);
title('Moderna');
subplot(2,2,4);
parcorr(PFE_pct);
title('Pfizer');

figure();
subplot(2,2,1);
qqplot(AZN_pct);
title('AstraZeneca');
subplot(2,2,2);
qqplot(JNJ_pct);
title('Johnson & Johnson');
subplot(2,2,3);
qqplot(MRNA_pct);
title('Moderna');
subplot(2,2,4);
qqplot(PFE_pct);
title('Pfizer');

% Plot first differences

clear g
g(1,1)=gramm('x',1:251,'y',AZN_pct);
g(1,1).geom_line();
g(1,1).set_names('x','Day','y','AZN');
g(1,1).axe_property('YLim',[-0.2 0.2]);
g(1,1).set_title('AstraZeneca');
g(1,2)=gramm('x',1:251,'y',JNJ_pct);
g(1,2).geom_line();
g(1,2).axe_property('YLim',[-0.2 0.2]);
g(1,2).set_names('x','Day','y','JNJ');
g(1,2).set_title('Johnson & Johnson');
g(2,1)=gramm('x',1:251,'y',MRNA_pct);
g(2,1).geom_line();
g(2,1).axe_property('YLim',[-0.2 0.2]);
g(2,1).set_names('x','Day','y','MRNA');
g(2,1).set_title('Moderna');
g(2,2)=gramm('x',1:251,'y',PFE_pct);
g(2,2).geom_line();
g(2,2).set_names('x','Day','y','PFE');
g(2,2).set_title('Pfizer');
g(2,2).axe_property('YLim',[-0.2 0.2]);
figure('Position',[100 100 800 800]);
g(1,1).set_color_options('map','d3_10');
g(1,2).set_color_options('map','brewer3');
g(2,1).set_color_options('map','winter');
g(2,2).set_color_options('map','brewer_pastel');
g.draw();

% histograms
clear g;
g2(1,1)=gramm('x',AZN_pct);
g2(1,2)=gramm('x',JNJ_pct);
g2(2,1)=gramm('x',MRNA_pct);
g2(2,2)=gramm('x',PFE_pct);
g2(1,1).stat_bin('normalization','probability','geom','overlaid_bar');
g2(1,2).stat_bin('normalization','probability','geom','overlaid_bar');
g2(2,1).stat_bin('normalization','probability','geom','overlaid_bar');
g2(2,2).stat_bin('normalization','probability','geom','overlaid_bar');
g2(1,1).set_title('AstraZeneca');
g2(1,2).set_title('Johnson & Johnson');
g2(2,1).set_title('Moderna');
g2(2,2).set_title('Pfizer');
g2(1,1).set_color_options('map','d3_10');
g2(1,2).set_color_options('map','brewer3');
g2(2,1).set_color_options('map','winter');
g2(2,2).set_color_options('map','brewer_pastel');
figure('Position',[100 100 800 600]);
g2.draw();

%% TEST STATIONARITY
% As before I test for unit root. Since all the p-values are low I can assume 
% stationarity of the time series.

%No TS

[~,pValue] = adftest(AZN_pct); 
pValue
[~,pValue] = adftest(JNJ_pct); 
pValue
[~,pValue] = adftest(MRNA_pct); 
pValue
[~,pValue] = adftest(PFE_pct); 
pValue


%TS
[~,pValue] = adftest(AZN_pct,'model','TS'); 
pValue
[~,pValue] = adftest(JNJ_pct,'model','TS'); 
pValue
[~,pValue] = adftest(MRNA_pct,'model','TS'); 
pValue
[~,pValue] = adftest(PFE_pct,'model','TS'); 
pValue

%% ARMA order for AZN
% We plot ACF and PACF and I observe a clear drop in both after the first lag. 
% So by exploratory analysis we would choose an ARMA(1,1). Then we proceed by 
% simulating the BIC and AIC coefficients. BIC identifies the best model in (1,1), 
% while AIC prefers the more complex (3,2) model. Overall we chose the (1,1) model.

LOGL = zeros(5,5); % Initialize
PQ = zeros(5,5);

 for p = 0:4
     for q = 0:4
        Mdl = arima(p,0,q);
        [EstMdl,~,logL] = estimate(Mdl,AZN_pct,'Display','off');
        LOGL(p+1,q+1) = logL;
        PQ(p+1,q+1) = p + q;
     end
end

LOGL = reshape(LOGL,25,1);
PQ = reshape(PQ,25,1);
[aic,bic] = aicbic(LOGL,PQ+1,100);
bicm=reshape(bic,5,5); %minimal value in (1,1)
aicm=reshape(aic,5,5); %minimal value in (3,3), but also (1,1) good value

%%

T=length(AZN_pct);
mod=arima(1,0,0);
AZN_mod=estimate(mod,AZN_pct);

AZN_res= infer(AZN_mod,AZN_pct);
AZN_pred = AZN_pct-AZN_res;

figure()
plot(AZN_pct,'-b')
hold on
plot(AZN_pred,'-r')

dy=AZN_pct;
Mdl = mod; % define the model
EstMdl = AZN_mod; % perform estimation
[dyF,dyMSE] = forecast(EstMdl,1,dy(end-10:end)); % forecast h step ahead %try last 10 days
dyF_ci = [dyF - 1.96*sqrt(dyMSE) dyF + 1.96*sqrt(dyMSE)]
figure
plot(dy,'Color',[.7,.7,.7]);
hold on
plot(T+1,dyF,'bo','LineWidth',2);
plot(T+1,dyF + 1.96*sqrt(dyMSE),'rd',...
		'LineWidth',2);
plot(T+1,dyF - 1.96*sqrt(dyMSE),'rd','LineWidth',2);
legend('Observed','Forecast',...
		'95\% Confidence Interval','Location','Best','interpreter','latex');
title('1-Step Ahead','interpreter','latex')
grid on
set(gca,'FontSize',20)

%% ARMA order for JNJ
% We plot ACF and PACF and I observe a clear drop in both after the first lag. 
% So by exploratory analysis we would choose an ARMA(1,1).  BIC and AIC both identify 
% the best model in (3,3).Overall we chose the (3,3) model.
LOGL = zeros(5,5); % Initialize
PQ = zeros(5,5);

 for p = 0:4
     for q = 0:4
        Mdl = arima(p,0,q);
        [EstMdl,~,logL] = estimate(Mdl,JNJ_pct,'Display','off');
        LOGL(p+1,q+1) = logL;
        PQ(p+1,q+1) = p + q;
     end
end

LOGL = reshape(LOGL,25,1);
PQ = reshape(PQ,25,1);
[aic,bic] = aicbic(LOGL,PQ+1,100);
bicm=reshape(bic,5,5); 
aicm=reshape(aic,5,5); 

%%
mod=arima(3,0,3);
JNJ_mod=estimate(mod,JNJ_pct);

JNJ_res= infer(JNJ_mod,JNJ_pct);
JNJ_pred = JNJ_pct-JNJ_res;
figure()
plot(JNJ_pct,'-b')
hold on
plot(JNJ_pred,'-r')

dy=JNJ_pct;
Mdl = mod; % define the model
EstMdl = JNJ_mod; % perform estimation
[dyF,dyMSE] = forecast(EstMdl,1,dy(end-10:end)); % forecast h step ahead 
dyF_ci = [dyF - 1.96*sqrt(dyMSE) dyF + 1.96*sqrt(dyMSE)]
figure
plot(dy,'Color',[.7,.7,.7]);
hold on
plot(T+1,dyF,'bo','LineWidth',2);
plot(T+1,dyF + 1.96*sqrt(dyMSE),'rd',...
		'LineWidth',2);
plot(T+1,dyF - 1.96*sqrt(dyMSE),'rd','LineWidth',2);
legend('Observed','Forecast',...
		'95\% Confidence Interval','Location','Best','interpreter','latex');
title('1-Step Ahead','interpreter','latex')
grid on
set(gca,'FontSize',20)

%% ARMA order for MRNA
% We plot ACF and PACF and we observe a quite clear behaviour after the first 
% lag. BIC and AIC both identify the best model in (2,2).Overall we chose the 
% (2,2).


LOGL = zeros(5,5); % Initialize
PQ = zeros(5,5);

 for p = 0:4
     for q = 0:4
        Mdl = arima(p,0,q);
        [EstMdl,~,logL] = estimate(Mdl,MRNA_pct,'Display','off');
        LOGL(p+1,q+1) = logL;
        PQ(p+1,q+1) = p + q;
     end
end

LOGL = reshape(LOGL,25,1);
PQ = reshape(PQ,25,1);
[aic,bic] = aicbic(LOGL,PQ+1,100);
bicm=reshape(bic,5,5); 
aicm=reshape(aic,5,5); 

%%
mod=arima(2,0,2);
MRNA_mod=estimate(mod,MRNA_pct); %low p-values for AR

MRNA_res= infer(MRNA_mod,MRNA_pct);
MRNA_pred = MRNA_pct-MRNA_res;
figure();
plot(MRNA_pct,'-b');
hold on;
plot(MRNA_pred,'-r');

dy=MRNA_pct;
Mdl = mod; % define the model
EstMdl = MRNA_mod; % perform estimation
[dyF,dyMSE] = forecast(EstMdl,1,dy(end-10:end)); % forecast h step ahead

dyF 
dyF_ci = [dyF - 1.96*sqrt(dyMSE) dyF + 1.96*sqrt(dyMSE)]

figure
plot(dy,'Color',[.7,.7,.7]);
hold on
plot(T+1,dyF,'bo','LineWidth',2);
plot(T+1,dyF + 1.96*sqrt(dyMSE),'rd',...
		'LineWidth',2);
plot(T+1,dyF - 1.96*sqrt(dyMSE),'rd','LineWidth',2);
legend('Observed','Forecast',...
		'95\% Confidence Interval','Location','Best','interpreter','latex');
title('1-Step Ahead','interpreter','latex')
grid on
set(gca,'FontSize',20)
%% ARMA order for PFE
% We plot ACF and PACF and I observe a not so clear behaviour after the first 
% lag. So by exploratory analysis we would choose an ARMA(1,1).
% BIC identifies the best model in (1,1), while AIC chooses (4,4).Overall we 
% chose the (1,1) model.

LOGL = zeros(4,4); % Initialize
LOGL = zeros(5,5); % Initialize
PQ = zeros(5,5);

 for p = 0:4
     for q = 0:4
        Mdl = arima(p,0,q);
        [EstMdl,~,logL] = estimate(Mdl,PFE_pct,'Display','off');
        LOGL(p+1,q+1) = logL;
        PQ(p+1,q+1) = p + q;
     end
end

LOGL = reshape(LOGL,25,1);
PQ = reshape(PQ,25,1);
[aic,bic] = aicbic(LOGL,PQ+1,100);
bicm=reshape(bic,5,5); 
aicm=reshape(aic,5,5); 

%%
mod=arima(4,0,4); %better choice than (1,1)
PFE_mod=estimate(mod,PFE_pct); %low p-values

PFE_res= infer(PFE_mod,PFE_pct);
PFE_pred = PFE_pct-PFE_res;

figure()
plot(PFE_pct,'-b')
hold on
plot(PFE_pred,'-r')

dy=PFE_pct;
Mdl = mod; % define the model
EstMdl = PFE_mod; % perform estimation
[dyF,dyMSE] = forecast(EstMdl,1,dy(end-10:end)); % forecast h step ahead

dyF 
dyF_ci = [dyF - 1.96*sqrt(dyMSE) dyF + 1.96*sqrt(dyMSE)]

figure
plot(dy,'Color',[.7,.7,.7]);
hold on
plot(T+1,dyF,'bo','LineWidth',2);
plot(T+1,dyF + 1.96*sqrt(dyMSE),'rd',...
		'LineWidth',2);
plot(T+1,dyF - 1.96*sqrt(dyMSE),'rd','LineWidth',2);
legend('Observed','Forecast',...
		'95\% Confidence Interval','Location','Best','interpreter','latex');
title('1-Step Ahead','interpreter','latex')
grid on
set(gca,'FontSize',20)
%% MEAN MODEL DIAGNOSTIC for AZN

AZN_stdr = AZN_res/sqrt(AZN_mod.Variance);
figure
subplot(2,2,1)
plot(AZN_stdr)
title('Standardized Residuals')
subplot(2,2,2)
histogram(AZN_stdr,10)
title('Standardized Residuals')
subplot(2,2,3)
autocorr(AZN_stdr)  
subplot(2,2,4)
parcorr(AZN_stdr) 
qqplot(AZN_stdr);
[h,p] = kstest(AZN_stdr);
p
[h,pValue] = lbqtest(AZN_res); 
pValue 
[h,pValue]=archtest(AZN_res);
pValue
% This tests the null hypothesis of no ARCH effects against the alternative ARCH model with k lags

%% MEAN MODEL DIAGNOSTIC for JNJ

JNJ_stdr = JNJ_res/sqrt(JNJ_mod.Variance);
figure
subplot(2,2,1)
plot(JNJ_stdr)
title('Standardized Residuals')
subplot(2,2,2)
histogram(JNJ_stdr,10) %skewed
title('Standardized Residuals')
subplot(2,2,3)
autocorr(JNJ_stdr)  
subplot(2,2,4)
parcorr(JNJ_stdr) 
qqplot(JNJ_stdr);
[h,p] = kstest(JNJ_stdr); 
p
[h,pValue] = lbqtest(JNJ_res); 
pValue
[h,pValue]=archtest(JNJ_res); 
pValue

%% MEAN MODEL DIAGNOSTIC for MRNA

MRNA_stdr = MRNA_res/sqrt(MRNA_mod.Variance);
figure
subplot(2,2,1)
plot(MRNA_stdr)
title('Standardized Residuals')
subplot(2,2,2)
histogram(MRNA_stdr,10)
title('Standardized Residuals')
subplot(2,2,3)
autocorr(MRNA_stdr)  
subplot(2,2,4)
parcorr(MRNA_stdr) 
qqplot(MRNA_stdr);
[h,p] = kstest(MRNA_stdr);
p
[h,pValue] = lbqtest(MRNA_res); 
pValue
[h,pValue]=archtest(MRNA_res); 
pValue
%% MEAN MODEL DIAGNOSTIC for PFE

PFE_stdr = PFE_res/sqrt(PFE_mod.Variance);
figure
subplot(2,2,1)
plot(PFE_stdr)
title('Standardized Residuals')
subplot(2,2,2)
histogram(PFE_stdr,10) 
title('Standardized Residuals')
subplot(2,2,3)
autocorr(PFE_stdr)  
subplot(2,2,4)
parcorr(PFE_stdr) 
qqplot(PFE_stdr);
[h,p] = kstest(PFE_stdr); 
p
[h,pValue] = lbqtest(PFE_res); 
pValue
[h,pValue]=archtest(PFE_res); 
pValue

%% PLOTS res sq

figure();
subplot(2,2,1);
autocorr(AZN_res.^2);
title('AstraZeneca');
subplot(2,2,2);
autocorr(JNJ_res.^2);
title('Johnson & Johnson');
subplot(2,2,3);
autocorr(MRNA_res.^2);
title('Moderna');
subplot(2,2,4);
autocorr(PFE_res.^2);
title('Pfizer');

figure();
subplot(2,2,1);
parcorr(AZN_res.^2);
title('AstraZeneca');
subplot(2,2,2);
parcorr(JNJ_res.^2);
title('Johnson & Johnson');
subplot(2,2,3);
parcorr(MRNA_res.^2);
title('Moderna');
subplot(2,2,4);
parcorr(PFE_res.^2);
title('Pfizer');


%% GARCH model for AZN
% We observed the ARCH effect test in 3 cases has very high p-value, however 
% we still chose to model the conditional volatility.


LOGL = zeros(5,4); 
PQ = zeros(5,4);
for p = 0:4
    for q = 1:4
        Mdl3 = arima(1,0,0); 
        CVarMdl3 = garch(p,q);
        Mdl3.Variance = CVarMdl3;
        [EstMdl,~,logL] = estimate(Mdl3,AZN_pct,'Display','off');
        LOGL(p+1,q) = logL;
        PQ(p+1,q) = p + q;
     end
end

LOGL = reshape(LOGL,20,1);
PQ = reshape(PQ,20,1);
[aic,bic] = aicbic(LOGL,PQ+1,100);
bicm=reshape(bic,5,4); 
aicm=reshape(aic,5,4);

%% 

dy=AZN_pct;
Mdl3 = arima(1,0,0); 
CVarMdl3 = garch(0,2);
Mdl3.Variance = CVarMdl3;
[GARCH_A,~,logL] = estimate(Mdl3,dy); 

%% GARCH model for JNJ

figure
subplot(1,2,1)
autocorr(JNJ_res.^2); 
subplot(1,2,2)
parcorr(JNJ_res.^2); 


LOGL = zeros(4,4); 
PQ = zeros(4,4);
for p = 1:4
    for q = 1:4
        Mdl3 = arima(3,0,3); 
        CVarMdl3 = garch(p,q);
        Mdl3.Variance = CVarMdl3;
        [EstMdl,~,logL] = estimate(Mdl3,JNJ_pct,'Display','off');
        LOGL(p,q) = logL;
        PQ(p,q) = p + q;
     end
end

LOGL = reshape(LOGL,16,1);
PQ = reshape(PQ,16,1);
[aic,bic] = aicbic(LOGL,PQ+1,100);
bicm=reshape(bic,4,4); 
aicm=reshape(aic,4,4); 
%% 


dy=JNJ_pct;
Mdl3 = arima(3,0,3); 
CVarMdl3 = garch(2,1);
Mdl3.Variance = CVarMdl3;
[GARCH_J,~,logL] = estimate(Mdl3,dy); 

%% GARCH model for MRNA

subplot(1,2,1)
autocorr(MRNA_res.^2);
subplot(1,2,2)
parcorr(MRNA_res.^2);
 

LOGL = zeros(5,4); 
PQ = zeros(5,4);
for p = 0:4
    for q = 1:4
        Mdl3 = arima(2,0,2); 
        CVarMdl3 = garch(p,q);
        Mdl3.Variance = CVarMdl3;
        [EstMdl,~,logL] = estimate(Mdl3,MRNA_pct,'Display','off');
        LOGL(p+1,q) = logL;
        PQ(p+1,q) = p + q;
     end
end

LOGL = reshape(LOGL,16+4,1);
PQ = reshape(PQ,16+4,1);
[aic,bic] = aicbic(LOGL,PQ+1,100);
bicm=reshape(bic,4+1,4); 
aicm=reshape(aic,4+1,4); 
%% 
dy=MRNA_pct;
Mdl3 = arima(2,0,2); 
CVarMdl3 = garch(0,2);
Mdl3.Variance = CVarMdl3;
GARCH_M = estimate(Mdl3,dy);


%% GARCH model for PFE

subplot(1,2,1)
autocorr(PFE_res.^2);
subplot(1,2,2)
parcorr(PFE_res.^2);

dimp=9;
dimq=9;
LOGL = zeros(dimp,dimq); 
PQ = zeros(dimp,dimq);
for p = 0:dimp-1
    for q = 1:dimq
        display(p);
        display(q);
        Mdl3 = arima(4,0,4); 
        CVarMdl3 = garch(p,q);
        Mdl3.Variance = CVarMdl3;
        [EstMdl,~,logL] = estimate(Mdl3,PFE_pct,'Display','off');
        LOGL(p+1,q) = logL;
        PQ(p+1,q) = p + q;
     end
end

LOGL = reshape(LOGL,dimp*dimq,1);
PQ = reshape(PQ,dimp*dimq,1);
[aic,bic] = aicbic(LOGL,PQ+1,100);
bicmx=reshape(bic,dimp,dimq); 
aicmx=reshape(aic,dimp,dimq);
%% 

dy=PFE_pct;
Mdl3 = arima(4,0,4); 
CVarMdl3 = garch(0,9); 
Mdl3.Variance = CVarMdl3;
GARCH_P = estimate(Mdl3,dy);

 
%% EGARCH

Mdl3 = arima(1,0,0); 
CVarMdl3 = egarch(1,1);
Mdl3.Variance = CVarMdl3;
[EGARCH_A,~,logL] = estimate(Mdl3,AZN_pct);


Mdl3 = arima(3,0,3); 
CVarMdl3 = egarch(1,1); %% bic 3,2 aic 4,2
Mdl3.Variance = CVarMdl3;
[EGARCH_J,~,logL] = estimate(Mdl3,JNJ_pct);


Mdl3 = arima(2,0,2); 
CVarMdl3 = egarch(1,1); %%best for aic and bic
Mdl3.Variance = CVarMdl3;
[EGARCH_M,~,logL] = estimate(Mdl3,MRNA_pct);


Mdl3 = arima(4,0,4); 
CVarMdl3 = egarch(1,1);
Mdl3.Variance = CVarMdl3;
[EGARCH_P,~,logL] = estimate(Mdl3,PFE_pct);


%% FORECASTING

%Load forecast data
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
AZNf = readtable('AZN_future.csv');
AZNf = AZNf(:,5);
AZNf = table2array(AZNf);
JNJf = readtable('JNJ_future.csv');
JNJf = JNJf(:,5);
JNJf = table2array(JNJf);
MRNAf = readtable('MRNA_future.csv');
MRNAf = MRNAf(:,5);
MRNAf = table2array(MRNAf);
PFEf = readtable('PFE_future.csv');
PFEf = PFEf(:,5);
PFEf = table2array(PFEf);


%Go to logs

AZNf=log(AZNf);
laggedAZNf = lagmatrix(AZNf,1);
diff=AZNf-laggedAZNf; 
AZNf_pct=diff;
AZNf_pct(1,:)=[];

JNJf=log(JNJf);
laggedJNJ = lagmatrix(JNJf,1);
diff=JNJf-laggedJNJ;
JNJf_pct=diff;
JNJf_pct(1,:)=[];

MRNAf=log(MRNAf);
laggedMRNAf = lagmatrix(MRNAf,1);
diff=MRNAf-laggedMRNAf;
MRNAf_pct=diff;
MRNAf_pct(1,:)=[];

PFEf=log(PFEf);
laggedPFEf = lagmatrix(PFEf,1);
diff=PFEf-laggedPFEf;
PFEf_pct=diff;
PFEf_pct(1,:)=[];

%Forecasting GARCH vs EGARCH

%A
[dy_AG,dyMSE_AG] = forecast(GARCH_A,10,AZNf_pct(1:10));
s=dy_AG-AZNf_pct(1:10);
MAE_AG=sum(abs(s))./10;
MSE_AG=sum(s.^2)./10;

[dy_AE,dyMSE_AE] = forecast(EGARCH_A,10,AZNf_pct(1:10));
s=dy_AE-AZNf_pct(1:10);
MAE_AE=sum(abs(s))./10;
MSE_AE=sum(s.^2)./10;

%J

[dy_JG,dyMSE_JG] = forecast(GARCH_J,10,JNJf_pct(1:10));
s=dy_JG-JNJf_pct(1:10);
MAE_JG=sum(abs(s))./10;
MSE_JG=sum(s.^2)./10;

[dy_JE,dyMSE_JE] = forecast(EGARCH_J,10,JNJf_pct(1:10));
s=dy_JE-JNJf_pct(1:10);
MAE_JE=sum(abs(s))./10;
MSE_JE=sum(s.^2)./10;

%M
[dy_MG,dyMSE_MG] = forecast(GARCH_M,10,MRNAf_pct(1:10));
s=dy_MG-MRNAf_pct(1:10);
MAE_MG=sum(abs(s))./10;
MSE_MG=sum(s.^2)./10;

[dy_ME,dyMSE_ME] = forecast(EGARCH_M,10,MRNAf_pct(1:10));
s=dy_ME-MRNAf_pct(1:10);
MAE_ME=sum(abs(s))./10;
MSE_ME=sum(s.^2)./10;

%P
[dy_PG,dyMSE_PG] = forecast(GARCH_P,10,PFEf_pct(1:10));
s=dy_PG-PFEf_pct(1:10);
MAE_PG=sum(abs(s))./10;
MSE_PG=sum(s.^2)./10;

[dy_PE,dyMSE_PE] = forecast(EGARCH_P,10,PFEf_pct(1:10));
s=dy_PE-PFEf_pct(1:10);
MAE_PE=sum(abs(s))./10;
MSE_PE=sum(s.^2)./10;



%% correlation
cmap = copper(6);

%subplot(3,2,1);

corr(AZN_pct,JNJ_pct);
[se1,CI] = bootstrap_R2(AZN_pct,JNJ_pct);
%h=hist(se1);
title('AZN_JNJ')
CI

%subplot(3,2,2);

corr(AZN_pct,MRNA_pct);
[se2,CI] = bootstrap_R2(AZN_pct,MRNA_pct);
%hist(se2);
title('AZN_MRNA')

CI

%subplot(3,2,3);

corr(AZN_pct,PFE_pct);
[se3,CI] = bootstrap_R2(AZN_pct,PFE_pct);
%hist(se3);
%title('AZN_PFE')

CI

%subplot(3,2,4);

corr(JNJ_pct,PFE_pct);
[se4,CI] = bootstrap_R2(JNJ_pct,PFE_pct);
%hist(se4);
%title('PFE_JNJ')

CI


%subplot(3,2,5);

corr(JNJ_pct,MRNA_pct);
[se5,CI] = bootstrap_R2(JNJ_pct,MRNA_pct);
%hist(se5);
CI
title('MRNA_JNJ')



%subplot(3,2,6);

corr(PFE_pct,MRNA_pct);
[se6,CI] = bootstrap_R2(PFE_pct,MRNA_pct);
%hist(se6);
title('PFE_MRNA')

CI

%%
figure();
g6(1,1)=gramm('x',se1);
g6(1,1).stat_bin('fill','face');
g6(1,1).set_title('''face''');
g6(1,1).set_title('AZN_JNJ');

g6(1,2)=gramm('x',se2);
g6(1,2).stat_bin('fill','face');
g6(1,2).set_title('''face''');
g6(1,2).set_title('AZN_MRNA');

g6(2,1)=gramm('x',se3);
g6(2,1).stat_bin('fill','face');
g6(2,1).set_title('''face''');
g6(2,1).set_title('AZN_PFE');

g6(2,2)=gramm('x',se4);
g6(2,2).stat_bin('fill','face');
g6(2,2).set_title('''face''');
g6(2,2).set_title('PFE_JNJ');

g6(3,1)=gramm('x',se5);
g6(3,1).stat_bin('fill','face');
g6(3,1).set_title('''face''');
g6(3,1).set_title('MRNA_JNJ');

g6(3,2)=gramm('x',se6);
g6(3,2).stat_bin('fill','face');
g6(3,2).set_title('''face''');
g6(3,2).set_title('PFE_MRNA');

g6(1,1).set_color_options('map','d3_10');
g6(1,2).set_color_options('map','brewer3');
g6(2,1).set_color_options('map','winter');
g6(2,2).set_color_options('map','brewer_pastel');
g6(3,1).set_color_options('map','brewer2');
g6(3,2).set_color_options('map','brewer1');


g6.draw();

