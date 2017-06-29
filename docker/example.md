# CEMiTool and Docker

##### The official Docker image for CEMiTool is publicly available on [Docker Hub](https://hub.docker.com/r/csblusp/cemitool/).
This example shows how to use _csblusp/cemitool_ Docker image with a complete R environment for running CEMiTool.
For more information about the CEMiTool analysis workflow, check the vignette.

### Initial consideration
All the examples here shown will consider a hypothetical `/home/me/myanalysis` directory which contains your _gene
expression_ and _sample annotation_ data. 
```bash
$ ls /home/me/myanalysis
expression.tsv  sample_annotation.tsv
```
Just to make things more practical, we will store this directory in a shell variable called `ANALYSIS_DIR`
```bash
$ export ANALYSIS_DIR=/home/me/myanalysis
```

### Running CEMiTool in interactive mode
With the following command you will be able to access the container in an interactive shell.
```bash
$ docker run --rm -it -v "$ANALYSIS_DIR:$ANALYSIS_DIR" -w $ANALYSIS_DIR csblusp/cemitool /bin/bash
root@1asd8asf1rf3:/home/me/myanalysis/$ 
```
Just to make sure you are in the right directory inside the container, take a look at the files in your current directory.

```bash
root@1asd8asf1rf3:/home/me/myanalysis/$ ls
expression.tsv  sample_annotation.tsv
```
Everything is just fine! Now you can start a R session, load the *CEMiTool* package and start your analysis.
```bash
root@1asd8asf1rf3:/home/me/myanalysis/$ R -q # Inside the container
> library(CEMiTool) # Load CEMiTool and perform your analysis
> expression_data <- read.table('expression.tsv', sep='\t')
> sample_annotation_data <- read.table('sample_annotation.tsv', sep='\t')
```
#### Output Files
As you are using [Docker volumes](https://docs.docker.com/engine/tutorials/dockervolumes/) (the **`-v`** flag), every file created
on `$ANALYSIS_DIR` directory inside the container will be available outside the container, in the same location. Let's first see
the files on our analysis directory
```bash
$ ls /home/me/myanalysis/
expression.tsv  sample_annotation.tsv
```
We just have two files, for now. Now we will create a new container and put a new file inside our `$ANALYSIS_DIR` inside the container

```bash
$ docker run --rm -it -v "$ANALYSIS_DIR:$ANALYSIS_DIR" -w $ANALYSIS_DIR csblusp/cemitool /bin/bash
root@4a2t3aasd3t4:/home/me/myanalysis/$ cat > output.txt
¯\_(ツ)_/¯
root@4a2t3aasd3t4:/home/me/myanalysis/$ exit
```
When you exit the container, everything is lost, except for the data persisted on volumes. Let's see if our _output.txt_ file 
created inside the killed container still exists.

```bash
$ ls $ANALYSIS_DIR
expression.tsv  sample_annotation.tsv  output.txt
$ cat $ANALYSIS_DIR/output.txt
¯\_(ツ)_/¯
```
Nice! You are now able to run your full co-expression analysis inside the _CEMiTool_ container.
