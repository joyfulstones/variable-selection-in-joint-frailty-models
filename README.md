# variable-selection-in-joint-frailty-models
This is the code and supplementary materials in the paper "Variable selection for joint frailty models of recurrent and terminal events".

We provide some SAS and Matlab code when the sample size n is 500 in this folder. Please change the paths in the code file before you run the program.

We only consider the case that the linear baseline hazard and log-Gamma frailties. The code for other settings are similar.

In the folder "SAS code", we provide the code for the full model method, the oracle method and the MIC method. For the full model method, the oracle method, and the MIC method, we only provide the code based on one replication (the data in the file "recurrent_event1.csv").  

Since the calculation by using bootstrap method takes too long, we conduct our simulation on the server. The code for bootstrap method are similar to that for the MIC method, which are not provided here.

In the subdirectory \code\matlab, we provide the code used to generate the data. You only need to run the program "main.m", then the original data are generated. In addition, bootstrap samples will be generated in the subdirectory\code\matlab by running the program "bs.m". "tempdata.m" is a function called by "main.m". Here we only provide the code that can generate one dataset and one bootstrap samples using one dataset.

In addition, we provide the SAS code that can be used to obtain the estimation of parameters using the oracle method. You only need to run the program "oracle.sas" in the subdirectory \code\sas, then the estimation of parameters are saved in a file named "estall.csv" in the subdirectory \code\sas code\result.

Similarly, we provide sas code used to obtain the estimation of parameters for the full model method and the MIC method. The operation is also under the subdirectory \code\sas code.


