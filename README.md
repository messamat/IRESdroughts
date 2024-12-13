# R code for 'Drought in intermittent river and ephemeral stream networks'

This repository contains R code associated with _Sarremejane, R., Messager, M. L., & Datry, T. (2022). 
Drought in intermittent river and ephemeral stream networks. Ecohydrology, 15(5), e2390. https://doi.org/10.1002/eco.2390_

## Abstract
Intermittent rivers and ephemeral streams (IRES), those watercourses that periodically
cease to flow or dry, are the world's most widespread type of river ecosystem.
Our understanding of the natural hydrology and ecology of IRES has greatly
improved, but their responses to extreme events such as drought remain a research
frontier. In this review, we present the state of the art, knowledge gaps and research
directions on droughts in IRES from an ecohydrological perspective. We clarify the
definition of droughts in IRES, giving recommendations to promote transferability in
how ecohydrological studies characterize droughts in non-perennial stream networks.
Based on a systematic search of the literature, we also identify common patterns
and sources of variation in the ecological responses of IRES to droughts and
provide a roadmap for further research to enable improved understanding and management
of IRES during those extreme hydrological events. Confusion in the terminology
and the lack of tools to assess the hydrological responses of IRES to drought
may have hindered the development of drought research in IRES. We found that
44% of studies confused the term drought with seasonal drying and that those that
measure droughts in a transferable way are a minority. Studies on ecological
responses to drought in IRES networks are still rare and limited to a few climatic
zones and organisms and mainly explored in perennial sections. Our review highlights
the need for additional research on this topic to inform IRES management and
conservation.

## Analysis structure and underlying data

**Structure**: the overall project directory is structured with the following sub-directories:  
data/ (raw data, read-only, not to be altered)  
results/ (results of the analysis, mostly reproduceable through code execution. However, also includes manually modified results)  
src/ (code written for the project)  

All scripts rely on this structure.

**Dependency management**: the R library of this project is managed by [renv](https://rstudio.github.io/renv/articles/renv.html).
This makes sure that the exact same package versions are used when recreating the project.
When calling `renv::restore()`, all required packages will be installed with their specific version. 
Please note that this project was built with R version 4.4 on a Windows 10 operating system.

**Syntax**: this analysis relies on the [data.table](https://rdatatable.gitlab.io/data.table/) syntax,
which provides a high-performance version of data.frame. 
It is concise, faster, and more memory efficient than conventional data.frames and the tidyverse syntax.

## Getting started
### Download the repository for R
In Git Bash, the following commands illustrate the procedure to make a local copy of the Github repository in a newly created directory at 
C://IRESdroughts/src :

```{r, engine = 'bash', eval = FALSE}
Mathis@DESKTOP MINGW64 /c/IRESdroughts/src
$ git clone https://github.com/messamat/IRESdroughts.git
```

In R Studio for Windows, the following procedure can be used:  

* Click on “File” in the menu ribbon  
* Select “New project…”  
* Choose the “Version control” option in the New Project Wizard window.
* Then, select “Git” in the next window.
* In the next window, fill the fields as follows:  
  * Repository URL: https://github.com/messamat/IRESdroughts
  * Project directory name: [will autofill as “IRESdroughts”]  
  * Create project as subdirectory of: [choose the parent directory of src]  
* Tick “Open in new session” and then click “Create project”.  

## Create the directory structure
In the IRESdroughts/ directory (parent of /src), create a /data and a /results directory.
Transfer /src/GRDCdat.zip in the data directory and unzip it.

## To reproduce figures  
To reproduce Figure 2, run plot_GRDC.R
