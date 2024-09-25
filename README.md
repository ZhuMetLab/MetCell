# MetCell

## Introduction
`MetCell` is an R package for end-to-end data processing for single-cell metabolomics. Now MetCell supports data processing from ion mobility-resolved mass cytometry (specifically timsTOF pro) with cell superposition and bottom-up assembly peak detection algorithm.

The docker image [`zhulab/metcell-r`](https://hub.docker.com/r/zhulab/metcell-r) contains entire environment for running `MetCell`. For convenience and taking fully use of `MetCell`, users can pull it and run `MetCell` just as following.

## What is metcell-r

`metcell-r` is a Docker environment to processing ion mobility-resolved mass cytometry data with `MetCell` R package. It is based on the [`r-base`](https://hub.docker.com/_/r-base/) docker.

## Pulling image

Users can pull the metcell-r image with the following script

```bash
docker pull zhulab/metcell-r
```

## Data preparation

The data folder to run MetCell should contain data files below: 

  * Raw data files (.d): Only one data file should be put in the file folder. And the file corresponding to the demo R script below could be download at the [deposit](https://www.biosino.org/node/sample/detail/OES00391980). (If you use this demo file, remember to UNZIP the data into .d format before you run MetCell)
  
  * R script: we provided demo [code](https://github.com/ZhuMetLab/MetCell/blob/main/extra/run.R). The R script must be named as "run.R". Parameters were described in the script. 
  
  * Time segment table: A .csv table recorded the time segment of raw file to be processed. The table must be named as "time_limit_table.csv".
  we provided demo [file](https://github.com/ZhuMetLab/MetCell/blob/main/extra/time_limit_table.csv). The colnames of the table must not be modified. The contain of 'data_file' column should be the same as your raw data file name. And start_time and end_time defined the time segment you want to process in the raw data, and the unit is second.
<br/><br>


## Run data processing work with metcell-r image

- go to your data folder (e.g., data)

```base
 cd data
```

- run docker using following code (*User should be permitted to run docker service*)

```bash
# MUST keep the code exactly as it is!
docker run -it --rm -v "$PWD":/data -u $(id -u ${USER}):$(id -g ${USER}) zhulab/metcell-r Rscript run.R
```

- wait till data processing work done

- `docker run` argument explanation:

    `-v "$PWD":/home/${USER}`: mapping current directory as home directory in docker container

    `-u $(id -u ${USER}):$(id -g ${USER})`: using current user to run the container

    `Rscript ~/run.R`: run run.R in container home directory with `Rscript`  command

## The result 

After the data processing work done, a folder name 'results' would be generated in the root folder. And the main results are listed following:

- "01_feature_table.csv": all detected peaks and their quantification in all individual cells;

- "02_isotope_annotation_table.csv": isotope annotation results for singly charged peaks and their quantification in all individual cells; 

- "03_metabolite_annotation_table.csv": putative metabolite annotation results for singly charged peaks through MS1 and CCS match with the library and their quantification in all individual cells.

# License
<a rel="license" href="https://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a>

This work is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)
