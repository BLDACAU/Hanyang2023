---
title: 'Exercise: Veteran Data'
author: "Jooyoung Lee"
date: "3/15/2023"
output:
  pdf_document:
    fig_caption: yes
    includes:
      in_header: header.tex
  word_document: default
  header-includes: \usepackage{booktabs}
  html_document:
    fig_caption: yes
  csl: biomed-central.csl
---


A study was conducted to compare the effects of two chemotherapy treatments in survival times for lung cancer patients. A total of 137 patients were randomly assigned to one of standard or test treatment group. The data include a number of covariates including tumor cell types, a Karnofsky performance score, the age of patient and time in months from diagnosis to randomization. 


\begin{table}[h!]
\caption{Description of Variables in Veteran Lung Cancer Data.}
\centering
\vspace{2mm}
\begin{tabular}{l l l}
\hline
Variable & Description\\
\hline
\textsf{trt} &  1 = standard 2 = test\\
\textsf{celltype}  &  1 = squamous, 2 = smallcell, 3 = adeno, 4 = large \\
\textsf{time} & survival time\\
\textsf{status} & event status\\
\textsf{karno} & Karnofsky performance score. A scale of 10.\\
\textsf{diagtime} & the months from diagnosis of lung cancer to entry into the study\\
\textsf{age} & the age of the patient in years\\
\textsf{prior} & 0 = no prior therapy 10 = prior therapy\\
\hline
\end{tabular}
\end{table}


```{r}
library(survival)
attach(veteran)
?veteran
head(veteran)
```

## 1. The Kaplan-Meier estimator

Now we calculate the Kaplan-Meier estimates of the survival function for each group using the \textsf{survfit} function in \textsf{library(survival)}. 

\textsf{survfit} uses the Greenwood's formula ($\tilde{\tau}^2(t)$ in Lecture note) for the variance calculation by default, which is 
$Var(\widehat{S}(t)) = \widehat{S}(t)^2 \sum\limits_{j:t_j \leq t} \frac{d_j}{n_j(n_j-d_j)}$, where $\widehat{S}(t) = \prod\limits_{j:t_j \leq t}\left(1-\frac{d_j}{n_j}\right)$.

### Confidence limits for $S(t)$
Here we consider $100(1-\alpha)\%$ pointwise confidence intervals for $S(t)$ for a particular specified time $t$. 

\textsf{survfit} has several options for types of confidence intervals. 
\begin{itemize}
\item[1.] \textsf{conf.type="plain"}
$$\widehat{S}(t) \pm z_{1-\alpha/2}\hspace{1mm} se(\widehat{S}(t))$$
\item[2.] \textsf{conf.type="log" (default)}
$$\widehat{S}(t) \exp(\pm z_{1-\alpha/2} \hspace{1mm} se(\log\widehat{S}(t)),$$
where 
\begin{align*}
Var(\log \widehat{S}(t)) =\sum_{j:t_j\leq t} \frac{d_j}{n_j(n_j-d_j)}.
\end{align*}

\item[3.] \textsf{conf.type="log-log"}
$$\{\widehat{S}(t)\}^{\exp[\pm z_{1-\alpha/2} \hspace{1mm} se(\log[-\log\{\widehat{S}(t)\}])]},$$
where
\begin{align*}
Var[\log \{-\log\widehat{S}(t)\}] &=\frac{1}{\{ \log \widehat{S}(t)\}^2} \sum_{j:t_j \leq t}\frac{d_j}{n_j(n_j-d_j)}.
 \end{align*}
\end{itemize}


```{r}
library(survival)
library(knitr)

# Kaplan-Meier curve for standard treatment 
fit0.km <- survfit(Surv(time, status) ~ 1, data=veteran, subset=(trt==1))  
summary(fit0.km)
fit0.km$surv   # est of S(t)
fit0.km$std.err  # se of H(t) or -log(S(t))

fit0.km$lower # lower 95% CI for S(t)
fit0.km$upper # upper 95% CI for S(t)

# Kaplan-Meier curve for test treatment 
fit1.km <- survfit(Surv(time, status) ~ 1, data=veteran, subset=(trt==2))  
summary(fit1.km)
fit1.km$surv   # est of S(t)
fit1.km$std.err  # se of H(t) or -log(S(t))

fit1.km$lower # lower 95% CI for S(t)
fit1.km$upper # upper 95% CI for S(t)
```


```{r, fit.height=2.5, fig.cap="Pointwise 95% interval estimates for S(t)"}
par(mfrow=c(1,2))
plot(fit0.km, xlab="TIME(DAYS) SINCE RANDOMIZATION", 
     ylab="ESTIMATED  PROBABILITY  OF  SURVIVAL", cex.lab=.7, cex.axis=0.5)
mtext("STANDARD TREATMENT", side=3, line=0.5,  cex=.8)
plot(fit1.km, xlab="TIME(DAYS) SINCE RANDOMIZATION", 
     ylab="ESTIMATED  PROBABILITY  OF  SURVIVAL", cex.lab=.7, cex.axis=0.5)
mtext("TEST TREATMENT", side=3, line=0.5,  cex=.8)
```

```{r, fit.height=3, fig.cap="Esimated survival function stratified by treatment group"}
library(ggsurvfit)
library(ggplot2)
library(ggsci)

veteran$trt2 <- factor(veteran$trt, 1:2, c("Standard", "Test"))

p1 <- survfit2(Surv(time, status) ~ trt2, data = veteran) %>%
  ggsurvfit(type="survival",
    theme = theme_classic() +
             theme(legend.position = "top") +
             theme(panel.grid.major.y = element_line(color = "gray90", size = 0.3))) +
  scale_x_continuous(breaks=c(0, 200, 400, 600, 800, 1000), expand=c(0,0))+
  scale_y_continuous(limits=c(0, 1), expand=c(0,0)) +
  labs(x = "TIME(DAYS) SINCE RANDOMIZATION", y = "SURVIVAL PROBABILITY",
       title = "MORTALITY") +
  theme(axis.title.y = element_text(vjust = -0.1),
        axis.title.x = element_text(vjust = -0.1),
        plot.title = element_text(hjust = -0.08, vjust=-3, size=10))+
  scale_color_jama()+
  add_risktable(times = c(0, 200, 400, 600, 800, 1000),
                risktable_stats = c("n.risk"),
                risktable_group = c("risktable_stats"),
                risktable_height = 0.2, # Adjusts the height of the risk table (default is 0.25),
                stats_label = "No. at risk",
                size =3, hjust= 0) 
   
p1
```

### Median Survival Time

Suppose we construct the confidence limit of $S(t)$ based on the log-transformation. From the graphical method, the median survival time for the standard treatment group is 103 days and the corresponding 95\% confidence intervals are (59, 132).  In addition, the median survival time for the test treatment group is 52 days and the corresponding 95\% confidence intervals are (44, 95).
```{r, fit.height=4, fig.align="center"}
par(mfrow=c(1,2))
plot(fit0.km, xlab="TIME(DAYS) SINCE RANDOMIZATION", 
     ylab="ESTIMATED  PROBABILITY  OF  SURVIVAL", cex.lab=0.7, cex.axis=0.5)
abline(h=0.5, lty=3, lwd=2, col="red")
mtext("STANDARD TREATMENT", side=3, line=0.5,  cex=.8) 

plot(fit1.km, xlab="TIME(DAYS) SINCE RANDOMIZATION", 
     ylab="ESTIMATED  PROBABILITY  OF  SURVIVAL", cex.lab=0.7, cex.axis=0.5)
abline(h=0.5, lty=3, lwd=2, col="red")
mtext("TEST TREATMENT", side=3, line=0.5,  cex=.8)

```


## 3. Log-rank test

Now, we want to test $H_0: S_1(t) = S_2(t)$ for all $t$.
```{r}
library(KMsurv)
survdiff(Surv(time, status)~factor(trt), data = veteran)

1-pchisq(0.00823, 1)
```
The test statistic is 
$$ \frac{(O-E)^2}{V} \sim \chi_{1}^2 \quad \text{under }H_0 $$
The p-values is $P(\chi_{1}^2 > 0.00823) = 0.9277156$. Therefore we do not reject the null hypothesis, therefore, the survival function is not different between the standard treatment group and the test treatment group. 


# Cox Proportional Hazard Regression Model

## Estimation

We fit a Cox regression model including the treatment indicator (\textsf{trt}), cell types (\text{celltype}), performance status (\textsf{karno}), the months from diagnosis of lung cancer to entry into the study (\textsf{diagtime}), prior therapy (\textsf{prior}), and age(\textsf{age}).


```{r}
library(survival)
attach(veteran)


veteran$celltypef. <- factor(veteran$celltype, levels=c("large", "squamous", 
                                                        "smallcell", "adeno"))
fit1 <- coxph(Surv(time, status)~ factor(trt) + celltypef. + karno + diagtime 
              + factor(prior) + age , data=veteran, method="breslow")
summary(fit1)
```
```{r}
library(knitr)
xname <- c("Treatment", "Cell Type", "", "", "performance status", 
           "Time from Diagnosis to entry", "Prior Therapy", "Age")
compare <- c("Test vs. Standard", "Squamous vs. Large", "Small vs. Large", 
             "Adeno vs. Large", "1 unit increase", "1 month increase", 
             "Yes vs. No", "1 year increase")
res <- data.frame(xname, compare, summary(fit1)$coef[,2],  
                  paste("(", round(exp(confint(fit1))[,1], 3), ", ",
                        round(exp(confint(fit1))[,2], 3), ")", sep=""), summary(fit1)$coef[,c(4,5)])
colnames(res) <- c("Covariates", "Comparison",  "RR", "95% CI", "z", "p-value")
rownames(res) <- NULL

kable(res,
digits = c(3, 3, 3, 4),
caption="The results of Cox regression modeling for lung cancer patient survival")
```


```{r}
# drop treatment
fit5 <- coxph(Surv(time, status)~  celltypef. + karno, data=veteran, method="breslow")
summary(fit5)
```

We drop non-significant variables by backward elimination, which leads to the model including only cell type and performance status variables with the form of 
\begin{align*}
h(t|\mathbf{z}_i) &= h_0(t)\exp(\mathbf{z}_i'\boldsymbol{\beta})\\
&=h_0(t)\exp(z_{i1}\beta_1 + z_{i2}\beta_2 + z_{i3}\beta_3 + z_{i4}\beta_4),
\end{align*}
where $\mathbf{z}_i = (I(\text{cell-type = Squamous}), I(\text{cell-type = Small cell}), I(\text{cell-type = Adeno cell}), \text{performance status})'$, $\boldsymbol{\beta} = (\beta_1, \beta_2, \beta_3, \beta_4)'$.


```{r}

xname <- c("Cell Type", "", "", "performance status")
compare <- c("Squamous vs. Large", "Small vs. Large", "Adeno vs. Large", "1 unit increase")
res <- data.frame(xname, compare, summary(fit5)$coef[,2], 
                  paste("(", round(exp(confint(fit5))[,1], 3), ", ", 
                  round(exp(confint(fit5))[,2], 3), ")", sep=""), summary(fit5)$coef[,c(4,5)])
colnames(res) <- c("Covariates", "Comparison",  "RR", "95% CI", "z", "p-value")
rownames(res) <- NULL

kable(res,
digits = c(3, 3, 3, 4),
caption="The results of the final model using a Cox regression for lung cancer patient survival")
```


Here we have $\hat{\beta}_1 = -0.325$, $s.e.(\hat{\beta}_1) = 0.276$ and the 95\% CI for $\beta_1$ is $-0.325 \pm 1.96\cdot 0.276 = (-0.867, 0.217)$. The estimated relative risk (RR) is $\exp(\hat{\beta}_1) = 0.722$ and the 95\% CI for RR is $\exp(-0.325 \pm 1.96\cdot 0.276) = (0.420, 1.243)$, which is interpreted as the relative risk of death for individuals with squamous cell type versus large cell type given a fixed value of performance status. 

Since the performance status increase by 10 units, we now estimate the effect of 10-unit change in performance status on the risk of death, controlling for cell type, which is $10\hat{\beta}_4$.

$10\hat{\beta}_4 = -0.309$, $s.e.(10\hat{\beta}_4) = 10*s.e.(\hat{\beta}_4) =0.052$ and the 95\% CI for $10\beta_4$ is $-0.309 \pm 1.96\cdot 0.052 = (-0.411, -0.208)$. Then, the estimated RR associated with a 10 unit increase in performance status is 0.734 and its 95\% CI is $\exp(-0.309 \pm 1.96\cdot 0.052) = (0.663, 0.813)$.  

We now plot the estimated cumulative hazard function and the survival function. In \textsf{R}, these estimators are calculated for a hypothesis person with the average covariates. Therefore, it is strongly recommended to use \textsf{newdata} to specify covariate values. 
For example, we are interested in estimating the cumulative hazard function and the survival function for the individuals with squamous cell type and performance status of 50. 

```{r, fig.cap="The estimated cumulative hazard function and the survival function for the individulas with cell type of squamous and performance status of 50."}

# The breslow estimates
H0 <- basehaz(fit5)

# The cumulative hazard and survival function for the individual with specified covariates
mydata <- with(veteran, data.frame(celltypef.="squamous", karno=50))
H <- survfit(fit5, newdata=mydata, type="aalen")

par(mfrow=c(1,2))
plot(H, fun="cumhaz", main="Estimated Cumulative Hazard")
plot(H, main="Estimated Survival Function")
```


## Model diagnostics

Here, we check the proportional hazard assumptions. 

```{r fig.width=4, fig.height=4, fig.cap="The Cox-Snell residual plots", fig.align="center"}
# Cox-Snell Residuals
coxsnell <- veteran$status - resid(fit5, type="martingale")
fit.coxsnell <- survfit(Surv(coxsnell, veteran$status)~1)

plot(log(fit.coxsnell$time), log(fit.coxsnell$cumhaz), ylab="log Cumulative Hazard of Cox-snell Residuals", xlab="log t")
abline(0, 1, lty=2, lwd=1.5)
```

Although the Cox-Snell residual plots in Figure 2 show no evidence for the misspecified model, we will conduct the formal test for the PH assumption. Note that the Cox-snell residual plots are not very informative. 

\textsf{cox.zph} function is used for checking the PH assumptions with the argument \textsf{transform} for a functional form of $g(t)$. \textsf{transform = "identity"} corresponds to the identity transform $g(t) = t$,  \textsf{transform = "log"} to $\log(t)$, \textsf{transform = "rank"} to the rank of the event times and \textsf{transform = "km"} to the Kaplan-Meier estimates $\hat{S}(t)$.  In the output of \textsf{cox.zph}, \textsf{chisq} gives the test statistics with \textsf{df}, degrees of freedom and \textsf{p} gives the p-value. The last row of \textsf{GLOBAL} gives the global test of proportional hazards over all p Covariates. 

```{r, fig.width=5, fig.height=4, fig.cap="The scaled Schoenfeld residual plot of the cell types in the model with the Kaplan-Meier transformation (i.e. $\\hat{\\beta}_l + r_{lj}^*(\\hat{\\boldsymbol{\\beta}})$ versus $t_j$, j=1,2.)"}
zph.id.fit5 <- cox.zph(fit5, transform  = 'identity') # g(t) = t
zph.log.fit5 <- cox.zph(fit5, transform  = 'log') # g(t) = log
zph.km.fit5 <- cox.zph(fit5, transform  = 'km') # g(t) = S(t) from the Kaplan-Meier - default

zph.id.fit5
zph.log.fit5
zph.km.fit5

plot(zph.id.fit5[1])
```
```{r, fig.width=5, fig.height=4, fig.cap="The scaled Schoenfeld residual plot of the performance status in the model with the Kaplan-Meier transformation (i.e. $\\hat{\\beta}_l + r_{lj}^*(\\hat{\\boldsymbol{\\beta}})$ versus $t_j$, l=1,2.)"}

plot(zph.id.fit5[2])
```

Based on the identity, log, and the Kaplan-Meier scale, there is some evidence for nonproportionality for the cell types and performance status, which are the significant predictor in the Cox model. 
In particular, the left penal of Figure 4 shows that the upward trend ends around 100 days for the scaled Schoenfeld residual plot of the performance status, $\hat{\beta}_l + r_{lj}^*(\hat{\boldsymbol{\beta}})$ versus $t_j$, which identifies 100 days as a point to allow the hazard ratio to change. 

## Stratification

To remove the problem of non-proportionality, first we use a cell-type variable as a stratifying variable. Then, a new model has the form of 
\begin{align*}
h_k(t|\mathbf{z}_{ki}) &= h_{k0}(t)\exp(z_{ki4}\beta_1), \quad k=1,2,3,4,
\end{align*}
which allows the cell-type specific baseline hazard functions. 

\textsf{strata} is used to specify a stratifying variable in \textsf{coxph}.

```{r, fig.width=5, fig.height=6, fig.cap="Survival functions versus time stratified by cell types for the individual with performance status of the mean of the performance status"}
fit5.str <- coxph(Surv(time, status) ~ karno + strata(celltypef.), data=veteran, method="breslow")
km <- survfit(fit5.str)
km
#summary(km) # gives the Kaplan-Meier estimates with 95% CI

library(survminer)
## Survival function versus time stratified by cell types
ggsurvplot(km, data = veteran, risk.table = TRUE, conf.int=FALSE, risk.table.height = 0.3) 
```

## Time-dependent Covariate

Second, we create a time-dependent covariate to allow the effect of performance status to depend on time through the interaction term with time. Here, we let $g(t) = I(t>100)$.

Then, our model has the form of 
\begin{equation}
h_k(t|\mathbf{z}_{ki}) = h_{k0}(t)\exp(z_{ki4}\beta_1 + z_{ki4}*I(t>100)\beta_2), \quad k=1,2,3,4
\end{equation}

The interpretation of regression coefficients is give in the following table.
\begin{table}[ht!]
\centering
\begin{tabular}{lllll}
\hline
Period & Cell type & performance status & Hazard & Relative Risk\\
\hline
$(0, 100]$ & $k$ & $10 + c$ & $h_{k0}(t)\exp((10+c)\beta_1)$ & $\exp(10\beta_1)$\\
 & $k$ &  $c$&  $h_{k0}(t)\exp(c\beta_1)$ &\\
 \\
$(100, \infty]$ &  $k$ & $10 + c$ & $h_{k0}(t)\exp((10+c)(\beta_1 + \beta_2) )$ & $\exp(10(\beta_1 + \beta_2))$\\
 & $k$ &   $c$&$ h_{k0}(t)\exp(c(\beta_1 + \beta_2))$ &\\
 \hline
\end{tabular}
\end{table}

```{r}
# Code to create a dataframe with a time-dependent covariate, It,
# where It = 0 if time <= 100 days; 1 if > 100 days months
timedepeff.f <- function(indata, cutpoint) {
       outdata <- NULL
       for (i in 1:nrow(indata)) {
         time   <- indata$time[i]
         status <- indata$status[i]
         if ( time <= cutpoint) {
           estart  <- 0
           estop   <- time
           estatus <- status
           It <- 0 
          } else{
           estart  <- c(0, cutpoint)
           estop   <- c(cutpoint,  time)
           estatus <- c(0, status)
           It      <- c(0, 1)
         }
        nlen <- length(estart)
        karno <- rep(indata$karno[i], nlen)
        id  <- rep(i, nlen)
        celltypef. <- rep(indata$celltypef[i], nlen)
        outdata <- rbind(outdata, data.frame(id, estart, estop, estatus, It, karno, celltypef.))
    }
    outdata <- data.frame(outdata, row.names=c(1:nrow(outdata)))
    return(outdata)
}


# adding a time-dependent covariate 
veteran.tdeff <- timedepeff.f(indata=veteran, cutpoint=100)
dimnames(veteran.tdeff)[[2]] <- c("id","estart","estop","estatus","timecut","karno", "celltypef.")

# adding interaction with karno*timecut and stratify by celltype
final.fit <- coxph(Surv(estart, estop, estatus) ~ karno + karno:factor(timecut) 
                + strata(celltypef.), data=veteran.tdeff, method="breslow")
summary(final.fit)

zph.id.final.fit <- cox.zph(final.fit, transform="identity")
zph.log.final.fit <- cox.zph(final.fit, transform="log")
zph.km.final.fit <- cox.zph(final.fit, transform="km")
zph.id.final.fit
zph.log.final.fit
zph.km.final.fit
```
```{r, fig.width=5, fig.height=4, fig.cap="The scaled Schoenfeld residual plot of the performance status in the model with the Kaplan-Meier transformation (i.e. $\\hat{\\beta}_1 + r_{1j}^*(\\hat{\\boldsymbol{\\beta}})$ versus $t_j$.)"}
plot(zph.id.final.fit[1])
```

The p-values from the marginal and global tests do not give evidence against the PH assumption and the plot of the scaled Schoenfeld residuals do not suggest any time trend in the coefficients. Therefore, we choose the model (1) as our final model.

```{r}
# Relative risk for a 10 unit increase for the first 100 days
rr1 <- exp(10*final.fit$coefficients[1])
rr.CI.1 <- exp(10*(final.fit$coefficients[1] + c(-1.96, 1.96)*sqrt(final.fit$var[1,1])))
rr1
rr.CI.1

# Relative risk for a 10 unit increase for after 100 days
est2<- final.fit$coefficients[1] + final.fit$coefficients[2]
sd2 <- sqrt(final.fit$var[1,1] + 2*final.fit$var[1,2] + final.fit$var[2,2])

rr2 <- exp(10*est2)
rr.CI.2 <- exp(10*(est2 + c(-1.96, 1.96)*sd2))
rr2
rr.CI.2
z <- est2/sd2 # z-value; H_0: beta1 + beta2 = 0
2*(1-pnorm(abs(z))) # p-value
```

For the first 100 days, the estimated RR associated with a 10 unit increase in performance status is 0.642 and its 95\% CI is $(0.567, 0.727)$ with p-value = 0.002. However, there is only 1.1\% reduction in the risk of death for a 10 unit increase in performance status after 100 days (RR=0.989, 95\% CI = $(0.776, 1.260)$; $p=0.929$. 
This suggests that performance status is more beneficial to reducing the risk of death of a lung cancer patient during the early days after the randomization.  
