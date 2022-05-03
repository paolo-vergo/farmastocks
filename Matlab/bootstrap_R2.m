function [se,CI] = bootstrap_R2(v1,v2)

%INPUT
% v1,v2 = vectors of equal length
%
%OUTPUT
%se=vector of correlations obtained via bootstrap
%CI= 95\% confidence interval for correlation

m=fitlm(v2,v1);
fitted=m.Fitted;
res=table2array(m.Residuals(:,1));
R2=m.Rsquared.Ordinary;

se = (bootstrp(10000,@(bootr)corr(fitted+bootr,v2),res));
m=mean(se);
std_sol=std(se);
CI=[m-1.96*std_sol m m+1.96*std_sol];

end
