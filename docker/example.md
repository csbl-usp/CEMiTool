# CEMiTool and Docker

##### The official Docker image for CEMiTool is publicly available on [Docker Hub](https://hub.docker.com/r/csblusp/cemitool/).
This example shows how to use the _csblusp/cemitool_ Docker image with a complete R environment for running *CEMiTool*.
For more information about the *CEMiTool* analysis workflow, please check the vignette. If you don't have Docker installed, you can follow the steps for your platform [here](https://docs.docker.com/engine/installation/)

### Initial consideration
All examples shown here will consider a hypothetical `/home/me/myanalysis` directory which contains your _gene
expression_ and _sample annotation_ data. 
```bash
$ ls /home/me/myanalysis
expression.tsv  sample_annotation.tsv
```
To make things more practical, we will store this directory in a shell variable called `ANALYSIS_DIR`
```bash
$ export ANALYSIS_DIR=/home/me/myanalysis
```

## Running CEMiTool in interactive mode
With the following command you will be able to access the container in an interactive shell.
```bash
$ docker run --rm -it --entrypoint=/bin/bash -v "$ANALYSIS_DIR:$ANALYSIS_DIR" -w $ANALYSIS_DIR -u $(id -u) csblusp/cemitool
docker@1asd8asf1rf3:/home/me/myanalysis/$ 
```
Just to make sure you are in the right directory inside the container, take a look at the files in your current directory.

```bash
docker@1asd8asf1rf3:/home/me/myanalysis/$ ls
expression.tsv  sample_annotation.tsv
```
Everything is just fine! Now you can start a R session, load the *CEMiTool* package and start your analysis.
```bash
docker@1asd8asf1rf3:/home/me/myanalysis/$ R -q # Inside the container
> library(CEMiTool) # Load CEMiTool and perform your analysis
> expression_data <- read.table('expression.tsv', sep='\t')
> sample_annotation_data <- read.table('sample_annotation.tsv', sep='\t')
```
---
## Running CEMiTool from command line
*CEMiTool* can also be run from the command line, with no need to enter the container. The following command will show the documentation for the command line version of *CEMiTool*
```bash
$ docker run --rm -v "$ANALYSIS_DIR:$ANALYSIS_DIR" -w "$ANALYSIS_DIR" -u $(id -u) csblusp/cemitool --help
```
The simplest way to run the *CEMiTool* cli is providing the **expression file** and **--output** parameters.
```bash
$ docker run --rm -v "$ANALYSIS_DIR:$ANALYSIS_DIR" -w "$ANALYSIS_DIR" -u $(id -u) csblusp/cemitool expression.tsv --output=output_directory
```
Yeah, I know, this command is lengthy. If you are running on Linux or MacOS, you can set an alias like this:
```bash
$ alias cemitool="docker run --rm  -v "$(pwd):$(pwd)" -w "$(pwd)" -u $(id -u) csblusp/cemitool"
```
Alternativelly, you can download [this file](cemitool) and put it into your `/usr/local/bin` directory.
```bash
$ wget https://raw.githubusercontent.com/csbl-usp/CEMiTool/master/docker/cemitool
$ chmod +x cemitool
$ mv cemitool /usr/local/bin
```
With this you can run *CEMiTool* container just by calling `cemitool` on the command line. As you can see, the command is a little different. We changed `$ANALYSIS_DIR` to `$(pwd)` just to make the command more general. That being said, in order to run the analysis using this alias you must go to your analysis directory:
```bash
$ cd $ANALYSIS_DIR
$ cemitool expression.tsv --output=output_directory
```
---
Have fun! :whale:
