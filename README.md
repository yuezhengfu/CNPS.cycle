# CNPS.cycle
This is an R package for element cycle analysis using metagenomic data.
![Graphical abstract](https://github.com/yuezhengfu/CNPS.cycle/assets/39332214/fa656480-d3f6-416b-8874-2d975ba44fcf)

# Citation
Zhengfu Yue, Jing Zhang, Rui Kou, Tianshun Liu,  Shoushan Sheng, Ye Tao, Liang Zeng, Zelong Zhao, Qinfen Li, Jing Zhang, Yukun Zou. CNPS.cycle: Streaming Shortgun Metagenomic Data Analysis for Soil Biogeochemical Cycles. Journal of Hazardous Materials (Submitted)

# Install
CNPS.cycle is available on Github, you can install it by:
```{r}
#Local installation (recommended, R ≥ 4.2.0)
1.Download "pkgs_local.zip" your local computer from https://github.com/yuezhengfu/CNPS.cycle/releases/download/V1.0.0/pkgs_local.zip
2.Unzip the `pkgs_local.zip` file and set it as the current working directory.
3.install.packages("CNPSInstaller_0.1.1.zip", repos = NULL)
4.library("CNPSInstaller")
5.CNPS_Installer(use_local = T ) #Install dependencies
6.CNPS_Install(use_local = T,pkg_type = "default" )

#Online installation (R ≥ 4.2.0)
library(devtools) 
install_github("yuezhengfu/CNPS.cycle")
```
# User Guide
---Demonstration using built-in dataset---

Download https://github.com/yuezhengfu/CNPS.cycle/releases/download/V1.0.0/SampleData_AutomatedExecutionScript.Rmd, a complete demonstration script. After ensuring that CNPS.cycle is correctly installed, open SampleData_AutomatedExecutionScript.Rmd in RStudio. Click "Run all" to automatically execute the script using the built-in dataset. All results will be stored in a directory named "Result," located in the same path as SampleData_AutomatedExecutionScript.Rmd.

---Processing your own data---

Preparing your data in the same format as the built-in dataset is crucial. Based on user feedback, we recommend not opening your data in Excel beforehand due to its row limit, which can result in data loss. Instead, in the "Load package and built-in data" section of SampleData_AutomatedExecutionScript.Rmd, read your datasets one by one. Once the data is successfully loaded, click "Run all" to execute the script.
 
# Workflow
![Figure 1](https://github.com/yuezhengfu/CNPS.cycle/assets/39332214/1cbe8b07-1e90-4a4c-973e-5ab89d34a2a9)

# Result layout
![1234](https://github.com/yuezhengfu/CNPS.cycle/assets/39332214/875f9ff2-978d-41fd-9b52-f5056e706ef5)
