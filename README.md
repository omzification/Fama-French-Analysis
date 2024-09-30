### Intro

In this analysis, I will explore portfolio performance using the methodology introduced by [Fama and French (1992)](https://en.wikipedia.org/wiki/FamaFrench_three-factor_model) and replicating some steps of the [Fama and Macbeth regression procedure](https://www.journals.uchicago.edu/doi/pdf/10.1086/260061) before performing a simple backtest. 

I downloaded portfolio and factor datasets from [French's Data Library](http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html). I used the '25 Portfolios Formed on Size and Book-to-Market (5 x 5)' and the 'Fama/French 3 Factors' datasets.

Notes:
- Monthly frequencies are used, for the period Jan 1960 to Dec 2022, included.
- The portfolios are not excess returns. The risk free interest rate is in the 3-factor file. This was used to compute excess returns on the portfolios.
- The risk free rate is a monthly rate (i.e. not annualised).


### Summary Statistics



I compute summary statistics (mean, variance, skewness, and kurtosis) for each portfolio's excess return using the dataset, where portfoliosData is an $`{n}`$x25 matrix. After calculating these statistics, I generate plots for each (e.g., plot(mean(portfoliosData)), plot(var(portfoliosData))) and label the x-axis to match the portfolios. Additionally, I have recorded the statistics for 'SMALL LoBM', 'SMALL HiBM', 'BIG LoBM', and 'BIG HiBM' below.

<div align="center">


||’SMALL LoBM’|	’SMALL HiBM’|	’BIG LoBM’|	’BIG HiBM’|
|:---:|:---:|:---:|:---:|:---:|
|**mean**|0.2904266|1.091956|0.5323569|0.6680848|
|**var**|63.93017|38.86011|22.26546|31.53922|
|**skewness**|-0.007667036|0.118535789|-0.249287352|-0.294635953|
|**kurtosis**	|4.859557	|8.056107	|4.325651	|5.063227|
</div>

<div align="center">
  
![image](https://github.com/user-attachments/assets/b0489f92-ff55-4b24-a36b-054c9e8312d8)
</div>



These results, which have a sample size of 756, were used to compute the t-statistic for the mean of each portfolio's excess returns. The t-statistic is calculated as the sample mean divided by the standard error (standard deviation divided by the square root of the sample size). I  then compare these t-statistics to the critical value for a 5% significance level. This will help identify which portfolio's excess returns have a mean that is statistically significant (i.e., different from zero) at the 5% level.


<div align="center">
  
||’SMALL LoBM’|	’SMALL HiBM’|	’BIG LoBM’|	’BIG HiBM’|
|:---:|:---:|:---:|:---:|:---:|
|**t-stat**|	0.9987214	|4.8163052	|3.1020436|	3.270898|
|significant? (Y/N)|	N	|Y|	Y|	Y|
</div>

### Time Series Analysis

For each portfolio ($`{ p = 1, 2, \ldots, 25 }`$), consider the regression:



$${R^{(p)} - r_f = \alpha^{(p)} + \beta_m^{(p)} (R_m - r_f) + \beta_s^{(p)} \text{SMB} + \beta_v^{(p)} \text{HML} + \varepsilon^{(p)} }$$

where:

- $`{ R^{(p)} }`$ is the $`{ n \times 1 }`$ vector of returns on the $`{ p }`$th portfolio.
- $`{ R_m - r_f }`$ (Mkt-RF) is the market excess return.
- SMB is the size factor (small minus big).
- HML is the value factor (high minus low).
- $`{ \alpha^{(p)}, \beta_m^{(p)}, \beta_s^{(p)}, \beta_v^{(p)} }`$ are the regression coefficients.
- $`{ \varepsilon^{(p)} }`$ is the error term.

All vectors are n × 1. These columns are labeled in the 3-factor file. By putting the model in the usual form

$${Y = Xb_0 + \varepsilon}$$

where:

$$
Y = \begin{pmatrix}
R_1^{(1)} - r_{f,1} & R_1^{(2)} - r_{f,1} & \dots & R_1^{(25)} - r_{f,1} \\
R_2^{(1)} - r_{f,2} & R_2^{(2)} - r_{f,2} & \dots & R_2^{(25)} - r_{f,2} \\
\vdots & \vdots & \ddots & \vdots \\
R_n^{(1)} - r_{f,n} & R_n^{(2)} - r_{f,n} & \dots & R_n^{(25)} - r_{f,n}
\end{pmatrix},
$$

$$ X = \begin{pmatrix}
1 & R_{m,1} - r_{f,1} & SMB_1 & HML_1 \\
1 & R_{m,2} - r_{f,2} & SMB_2 & HML_2 \\
\vdots & \vdots & \vdots & \vdots \\
1 & R_{m,n} - r_{f,n} & SMB_n & HML_n
\end{pmatrix}, $$

$$
b_0 = \begin{pmatrix}
\alpha^{(1)} & \alpha^{(2)} & \dots & \alpha^{(25)} \\
\beta_m^{(1)} & \beta_m^{(2)} & \dots & \beta_m^{(25)} \\
\beta_s^{(1)} & \beta_s^{(2)} & \dots & \beta_s^{(25)} \\
\beta_v^{(1)} & \beta_v^{(2)} & \dots & \beta_v^{(25)}
\end{pmatrix},
$$

$$
\varepsilon = \begin{pmatrix}
\varepsilon_1^{(1)} & \varepsilon_1^{(2)} & \dots & \varepsilon_1^{(25)} \\
\varepsilon_2^{(1)} & \varepsilon_2^{(2)} & \dots & \varepsilon_2^{(25)} \\
\vdots & \vdots & \ddots & \vdots \\
\varepsilon_n^{(1)} & \varepsilon_n^{(2)} & \dots & \varepsilon_n^{(25)}
\end{pmatrix},
$$

the OLS estimator is saved to a 4 × 25 matrix, $`{B}`$.

Each column of the matrix $`{B′}`$ is plotted below separately.

<div align="center">
  
![image](https://github.com/user-attachments/assets/185689d2-5613-42b4-b138-097908efed48)
</div>

The residuals for each individual portfolio regression ($`{ p = 1, 2, \ldots, 25 }`$):

$$ \varepsilon^{(p)}:= Y^{(p)} - XB^{(p)} $$

where $`Y^{(p)}`$ is the $`{p^{th}}`$ column in $`Y`$ , and $`B^{(p)}`$ is the $`{p^{th}}`$ column in $`B`$
and an estimator for $`VarB^{(p)}`$􏰁were computed. Recorded below are the values of $`B^{(p)}`$ (the estimated alphas and betas) and their standard errors for the portfolios
corresponding to ’SMALL LoBM’, ’SMALL HiBM’, ’BIG LoBM’, ’BIG HiBM’.



<div align="center">

||’SMALL LoBM’	|’SMALL HiBM’	|’BIG LoBM’|	’BIG HiBM’|
|:---:|:---:|:---:|:---:|:---:|
|**$`{\displaystyle {\hat {\alpha }}}`$ value (s.e.)**|	-0.4539688 (0.09842763)|	0.1752701 (0.07040636)|	0.1427329 (0.04040566)|	-0.1944998 (0.09072885)|
|**$`\hat{\beta}_m`$ value (s.e.)**| 	1.0942033 (0.02294111)|	0.9447934 (0.01641002)|	0.989888 (0.009417585)|	1.1246657 (0.02114671)|
|**$`\hat{\beta}_s`$ value (s.e.)** |	1.4010419 (0.03410647)|	1.106314 (0.02439673)|	-0.2377553 (0.014001095)	|-0.1442474 (0.03143874)|
|**$`\hat{\beta}_v`$ value (s.e.)**|	-0.2520452 (0.03413614)|	0.6940601 (0.02441796)|	-0.3510167 (0.014013275)	|0.8537941 (0.03146609)|
</div>


Then a t-test is used to test whether the intercepts and the betas are significantly different from zero at the 5% significance level.


<div align="center">
  
||’SMALL LoBM’	|’SMALL HiBM’|	’BIG LoBM’|	’BIG HiBM’|
|:---:|:---:|:---:|:---:|:---:|
|**$`\hat{\alpha}`$ t-stat**	|-4.612209	|2.489407|	3.532497|	-2.143748|
|significant? (Y/N)| 	Y|	Y	|Y	|Y|
|**$`\hat{\beta}_m`$ t-stat**| 	47.696187|	57.574163|	105.110604	|53.183968|
|significant? (Y/N)	|Y|	Y	|Y	|Y|
|**$`\hat{\beta}_s`$ t-stat** |	41.078473|	45.346809|	-16.981191|	-4.588206|
|significant? (Y/N)|	Y	|Y	|Y	|Y|
|**$`\hat{\beta}_v`$ t-stat**| 	-7.383528	|28.42417	|-25.04887|	27.133781|
|significant? (Y/N)|	Y	|Y	|Y	|Y|
</div>

### Cross-Sectional Regression

Suppose that there are three scalars, such that, for each portfolio ($`{ p = 1, 2, \ldots, 25 }`$)

$$ {\mathbb{E}(R_t^{(p)}-r_f) = \beta_m^{(p)}\lambda_m + \beta_s^{(p)}\lambda_s + \beta_v^{(p)}\lambda_v }$$

The $`{\lambda}`$’s are the risk premia associated to the three (pricing) factors: market, size (SMB) and value (HML). To estimate them we can run the following regression. At time $`{t}`$ (the $`{t^{th}}`$ row in the dataset)  

$${ (R_t^{(p)}-r_f) = \alpha_t + \beta_m^{(p)}\lambda_{m,t} + \beta_s^{(p)}\lambda_{s,t} + \beta_v^{(p)}\lambda_{v,t} + Z_t^{(p)} }$$

Hence, supposing that this equation holds for each $`{p}`$ where $`{\alpha_t}`$ is zero when we average across $`{t}`$ (if (1) holds) . The term $`{Z_t^{(p)}}`$ is a mean zero pricing error. We want to estimate the unknown parameters $`{\lambda_{m,t}, \lambda_{s,t}, \lambda_{v,t}}`$ and $`{\alpha_t}`$. We do not know $`{\beta_m^{(p)}, \beta_s^{(p)}, \beta_v^{(p)}}`$, but we have estimated and saved them earlier in the matrix $`{B}`$. Suppose that the estimated betas in $`{B'}`$ are the matrix of explanatory variables and regress at time $`{t}`$, the 25 portfolio returns on columns 2 to 4 in $`{B'}`$; this is made clear next. To do so, let $`{R_t}`$ be the 25×1 vector of returns at time $`{t}`$ where the $`{p^{th}}`$ row corresponds to the $`{p^{th}}`$ portfolio. At each time $`{t}`$, run the regression supposing the true model

$${Y = X b_0 + \varepsilon,}$$

where:

$${Y = \begin{pmatrix}
R_t^{(1)} - r_{f,t} \\
R_t^{(2)} - r_{f,t} \\
\vdots \\
R_t^{(25)} - r_{f,t}
\end{pmatrix},}$$

$${X = \begin{pmatrix}
1 & \beta_m^{(1)} & \beta_s^{(1)} & \beta_v^{(1)} \\
1 & \beta_m^{(2)} & \beta_s^{(2)} & \beta_v^{(2)} \\
\vdots & \vdots & \vdots & \vdots \\
1 & \beta_m^{(25)} & \beta_s^{(25)} & \beta_v^{(25)}
\end{pmatrix},}$$

$${b_0 = \begin{pmatrix}
\alpha_t \\
\lambda_{m,t} \\
\lambda_{s,t} \\
\lambda_{v,t}
\end{pmatrix},}$$

$${\varepsilon = \begin{pmatrix}
\varepsilon_t^{(1)} \\
\varepsilon_t^{(2)} \\
\vdots \\
\varepsilon_t^{(25)}
\end{pmatrix}.}$$

The goal is to find the estimator $`{ b = \begin{pmatrix} \hat{\alpha}_t \\ \hat{\lambda}_{m,t} \\ \hat{\lambda}_{s,t} \\ \hat{\lambda}_{v,t} \end{pmatrix} }`$ for $`{ b_0 }`$. And to repeat the procedure for each $`{t}`$ and save $`\hat{b}`$ as the $`{t}`$ column in the 4 × n matrix riskPremia. Plotted below are the time series of risk premia for all three factors separately and the pricing errors $`{\alpha_t}`$'s. 

<div align="center">
  
![image](https://github.com/user-attachments/assets/cc8449e5-afa0-42dd-9600-5bb3794579de)
</div>

Then the following statistics for the three risk premia and the pricing error were computed.

<div align="center">
  
||$`\hat{\alpha}	`$|$`{\hat{\lambda}}_m `$|$`	{\hat{\lambda}}_s	`$|$`{\hat{\lambda}}_v`$|
|:---:|:---:|:---:|:---:|:---:|
|**mean**	|1.1678016|	-0.6002972|	0.1232887	|0.3403902|
|**var**|	49.746484	|69.570288|	9.488296	|8.990091|
|**skewness**|	0.1292128|	-0.1396704|	0.8604242	|0.0439532|
|**kurtosis**	|4.645125	|5.114728|	10.125036|	6.722539|
</div>

A t-statistic for the pricing error $`{\alpha_t}`$ was computed to see if the statistic were asymptotically normally distributed, would you see evidence of omitted factors at the 5% level of significance.

<div align="center">
  
||$`\hat{\alpha}`$|
|:---:|:---:|
|**t-stat**|4.552488|
|omitted factors? (Y/N)|	Y|
</div>

### Backtest

Using the excess returns on the 25 portfolios, the sample covariance matrix was computed on the sample period Jan 2010 to Dec 2017, denoted by $`{COV}`$. The shrunk covariance matrix was computed by $`{COVS=(1-a)*average(diag(COV))*Identity(25)+a*COV}`$ 
with $`{a=0.5}`$. 
By “$`{average(diag(COV))*Identity(25)}`$” I mean you average the diagonal entries in $`{COV}`$ and then multiply the 25 × 25 identity matrix by this value.

I then calculate the portfolio weights that solve the minimum variance portfolio problem, with weights summing to one, using the two different estimators of the covariance matrix and denote them by $`{w}`$ and $`{ws}`$. I also consider the equally weighted portfolio, and denote its weights by $`{wEW}`$. I carry out a backtest of the three portfolios. To do so, I must use data unseen by the estimators $`{w}`$ and $`{ws}`$. Hence, I used the testing period Jan 2018 to Dec 2022 to compute the portfolio excess returns for the strategy that uses w, ws and wEW as portfolio weights. (Conducting a backtest on data used for estimation is said to be subject to forward looking bias and must be avoided when possible.) Also, I used the excess returns on the market $`{(R_{m,t}-r_f)}`$ as a benchmark. I compute the mean, variance, skewness and kurtosis for the three portfolios and the benchmark.

<div align="center">
  
||	$`{w}`$|$`{ws}`$|	$`{wEW}`$|	benchmark|
|:---:|:---:|:---:|:---:|:---:|
|**mean**|	0.45928655	|0.7147997|	0.6898007	|0.7516667|
|**var**|	24.47927069	|30.1188906	|43.6113149	|31.383265|
|**skewness**	|0.09288124	|-0.4888421	|-0.437475|	-0.2955828|
|**kurtosis**|	3.43976676|	3.6301979	|4.5988648	|2.9437211|
</div>

I also plot the cumulative (monthly compounded) excess returns of the portfolios and the benchmark. The excess returns are in percentage points. Hence, I make sure that I divide them by 100 in order to compute the cumulative compounded returns so that for example, 0.15 means 15%. If $`{R_s}`$ is the one month simple return at time $`{s}`$, the cumulative compounded return at time $`{t}`$ is $`{C_t = \prod_{s = 1}^{t} (1+R_s)}`$. The plot is the graph of $`{C_t}`$ against $`{t}`$.

<div align="center">
  
![image](https://github.com/user-attachments/assets/b5373329-3c98-4754-bac6-d0b4726fd6fe)
</div>

The Sharpe ratio for a monthly investment is defined as:

$$\frac{\sqrt{12 }\hat{\mu}}{\sqrt{\hat{\sigma}^2}}$$ 

where $`{\hat{\mu}}`$ is the sample mean of the monthly returns, and $`{\hat{\sigma}^2}`$ is the sample variance of the monthly returns. I compute the Sharpe ratios for the four portfolio strategies and test whether they are significantly greater than zero at the 5% and 10% level. I assume that the returns are independent identically distributed with finite second moment. The number of observations in the backtest is 60.

<div align="center">
  
||	$`{w}`$|$`{ws}`$|	$`{wEW}`$|	benchmark|
|:---:|:---:|:---:|:---:|:---:|
|**Sharpe Ratio**	|0.3215697	|0.4511859	|0.3618384|	0.4648008|
|Significant at 10%? (Y/N)|	N	|N	|N	|N|
|Significant at 5%? (Y/N)|	N	|N|	N|	N|
</div>





