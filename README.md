# Benchmarking machine learning and Cox regression approaches for lung cancer prognosis on combined RNA-seq and personalised metabolic modelling data
Lung cancer rates are rising, with 29.5M cases and 16.4M deaths projected by 2040; in 2020, there were 19.3M cases, 2.21M related to lung cancer, motivating research into personalized treatment using multi-omics analysis, genome-scale metabolic models, and machine learning, identifying 17 new biomarkers and 4 pathways.

## Aims and Objectives
(i) To evaluate and compare the accuracy of traditional machine learning and semi-parametric survival analysis methods in predicting patient outcomes. The objective of this aim to collect the patient survival data from the relevant repository and compare the predictive accuracy of the two methods using appropriate statistical metrics. Also, to provide recommendations for the preferred machine learning survival method for anticipating patient outcomes based on comparative analysis.

(ii) To investigate the utility of a consolidated transcriptomic and fluxomic model for accurately forecasting overall patient survival and determine if it outperforms the use of either data type alone. The intention behind this aim to compare the performance of the consolidated model with that of the individual models and to determine whether the consolidated model provides more accurate forecasts of overall patient survival compared to using either data type alone.

(iii) To identify novel biomarkers and pathways associated with lung cancer adenocarcinoma and explore their potential applications in advancing personalised medicine for this patient outcomes. The purpose of this aim to perform omics data analysis (transcriptomic and fluxomic) on relevant patient samples to identify potential novel biomarkers and to conduct pathways analysis to identify significant pathways associated with lung cancer adenocarcinoma. Also, to validate the novel biomarkers and pathways findings from the literature review. Furthermore, to assess the potential clinical relevance and applicability of these novel biomarkers and pathways in personalised medicine for lung cancer adenocarcinoma patients.

## Steps to Run the DT, GBT, RF, SVR, Cox Regression Models Using Feature Selection Methods PCA and RF
(i) First, create a folder and provide a name like "Lungcancer."

(ii) Download the datasets (Flux.xlsx and Geneexpression.xlsx) from https://github.com/pradipms/Lungcancer/tree/master/MACHINE%20LEARNING/PYTHON/MODELLING folder. Save the datasets in the "Lungcancer/Dataset" folder.

(iii) Download the machine learning modeling codes for PCA and RF feature selection at https://github.com/pradipms/Lungcancer/tree/master/MACHINE%20LEARNING/PYTHON/MODELLING/PCA and https://github.com/pradipms/Lungcancer/tree/master/MACHINE%20LEARNING/PYTHON/MODELLING/RF. Save them in the "Lungcancer/PCA" and "Lungcancer/RF" folders.

(iv) The machine learning modeling codes are in .ipynb format. Install Jupyter Notebook and upload code files in Jupyter Notebook.

(v) First, set the working directory (provide the correct directory path in the code), and then run the code.

## Steps to Run R Code
Some of the plot analyses were done in R and R Studio software.

(i) First, create a folder "R" in the "Lungcancer" folder.

(ii) Download the datasets from https://github.com/pradipms/Lungcancer/tree/master/MACHINE%20LEARNING/R/Dataset folder. Save the datasets in the "Lungcancer/R/Dataset" folder.

(iii) Codes are provided at https://github.com/pradipms/Lungcancer/tree/master/MACHINE%20LEARNING/R folder. Save the code in the "Lungcancer/R" folder.

(iv) Download R Software with R Studio and install it.

(v) Open R files in R Studio, and before running the files, set the working directory.

## Workflow
![MAIN_DESIGN](https://github.com/pradipms/Lungcancer/assets/78509712/201bb0ee-0078-4db4-9da0-5d9ecf1642a9)
