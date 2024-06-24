# CNPS.cycle
This is an R package for element cycle analysis using metagenomic data.
![GF](https://github.com/yuezhengfu/CNPS.cycle/assets/39332214/183e531f-31ff-4bb0-9504-0635b67422a7)
# Citation
Zhengfu Yue, Jing Zhang, Rui Kou, Tianshun Liu, Ye Tao, Liang Zeng, Zelong Zhao, Shoushan Sheng, Qinfen Li, Jing Zhang, Yukun Zou.CNPS.cycle: An R package for element cycle gene and functional microbial analysis based on shotgun metagenomic data. Soil Biology & Biochemistry (Submitted)
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
# Workflow
![Figure 1](https://github.com/yuezhengfu/CNPS.cycle/assets/39332214/1cbe8b07-1e90-4a4c-973e-5ab89d34a2a9)
# Result layout
![1234](https://github.com/yuezhengfu/CNPS.cycle/assets/39332214/875f9ff2-978d-41fd-9b52-f5056e706ef5)
