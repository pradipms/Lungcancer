install.packages("readr")
library("readr")
# Load libraries
library("survival")
library("survminer")
library("devtools")
library("readr")

flux_data<- read_csv("C:/Users/U0033207PC/Desktop/GENE Enrichment/FINAL VISUALISATION/ANALYTICS/BIOMARKER_RISK/BIOMARKER_RISK.csv")

#Define variables
dt1 <- c('Hugo_Symbol', 'FACOAL161','FACOAL226','FAOXC180x','FAOXC200180x','LPS2e','FAOXC22C20x','FAOXC5C5OHm','FOLR2','MTHFC','MTHFD2','MTHFDm','MTHFR3','PGI','GAPD','LDH_L','PGK','PGM','Risk','Age', 'Cancer_Tumor_Stage', 'Sex', 'Patient\'s Vital Status','Overall Survival (Months)')
flux_data1 <- flux_data[, !names(flux_data) %in% dt1 ]


flux_data1 <- flux_data[ , apply(flux_data, MARGIN = 2, function(x) sum(is.na(x)) == 0)]

#Feature Scaling
flux_data1[,2:18] = scale(flux_data1[,2:18])   

# Compute survival curves For Risk
fit <- survfit(Surv(flux_data1$'Overall Survival (Months)', flux_data1$'Patient\'s Vital Status') ~ Risk, data = flux_data1)
print(fit)

# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table

#Access to the value returned by survfit()
d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower
)
head(d)

# Plot with Censor of Risk
low_high_risk_plot <- ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Overall Survival (Months)",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 200.
  ggtheme = theme_bw() + theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ) +  theme(axis.line = element_line(color = "black")), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = c("Low Risk", "High Risk"),    # change legend labels.
  palette = c("#E7B800", "#2E9FDF") # custom color palettes.
)

low_high_risk_plot 

# Save Risk Plot in Pdf
pdf("low_high_risk_plot.pdf")
print(low_high_risk_plot, newpage = FALSE)
dev.off()

# Compute survival curves For Sex
fit <- survfit(Surv(flux_data1$'Overall Survival (Months)', flux_data1$'Patient\'s Vital Status') ~ Sex, data = flux_data1)
print(fit)

# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table

#Access to the value returned by survfit()
d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower
)
head(d)

# Plot with Censor of Sex
sex_plot<-ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Overall Survival (Months)",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 200.
  ggtheme = theme_bw() + theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
  ) +  theme(axis.line = element_line(color = "black")), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Male", "Female"),    # change legend labels.
  palette = 
    c("#E700B8", "#E7B800") # custom color palettes.
)

sex_plot

# Save Sex Plot in Pdf
pdf("sex_plot.pdf")
print(sex_plot, newpage = FALSE)
dev.off()

# Compute survival curves For Cancer Tumor Stage
fit <- survfit(Surv(flux_data1$'Overall Survival (Months)', flux_data1$'Patient\'s Vital Status') ~ Cancer_Tumor_Stage, data = flux_data1)
print(fit)

# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table

#Access to the value returned by survfit()
d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower
)
head(d)

# Plot Cancer Tumor Stage
cancer_tumor_stage_plot<-ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           xlab = "Overall Survival (Months)",   # customize X axis label.
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() + theme(
             plot.background = element_blank()
             ,panel.grid.major = element_blank()
             ,panel.grid.minor = element_blank()
           ) +  theme(axis.line = element_line(color = "black")), # customize plot and risk table with a theme.
           palette = c("#FF337A", "#D433FF", "#8033FF", "#33F3FF", "#3393FF", "#33FFA8", "#33FF36", "#FFC133", "#FF6833"))

cancer_tumor_stage_plot

# Save Cancer Tumor Stage Plot in Pdf
pdf("cancer_tumor_stage_plot.pdf")
print(cancer_tumor_stage_plot, newpage = FALSE)
dev.off()

# Plot of 17 Biomarkers
biomarkers_plot<-ggsurvplot(
  fit = survfit(Surv(flux_data1$'Overall Survival (Months)', flux_data1$'Patient\'s Vital Status') ~ 1, data = flux_data1), 
  risk.table = TRUE, # Add risk table
  xlab = "Overall Survival (Months)",
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_bw() + theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
  ) +  theme(axis.line = element_line(color = "black")), # customize plot and risk table with a theme.
  palette = "Red"
 )

biomarkers_plot

# Save 17 Biomarkers Plot in Pdf
pdf("biomarkers_plot.pdf")
print(biomarkers_plot, newpage = FALSE)
dev.off()

# Plot with Censor of Risk
test_plot <- ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Overall Survival (Months)",   # customize X axis label.
  break.time.by = 50,     # break X axis in time intervals by 200.
  ggtheme = theme_bw() + theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
  ) +  theme(axis.line = element_line(color = "black")), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = c("All"),    # change legend labels.
  palette = c("Red") # custom color palettes.
)


