# MedZIM

## Installation

The package could be installed by following line: 
```
devtools::install_github("quranwu/MedZIM")
```

## Usage

The main function is `MedZIM_func()`. The function could be run with sample dataset as: 
```
res<-MedZIM_func(dat=data_ZIM, xVar = "x",yVar = "y1",taxon_name = "taxon",libSize_name = "libSize",paraJobs=2)
```

The results of each taxon could be accessed by: 
```
res$fullList$taxon1
res$fullList$taxon2
...
```
