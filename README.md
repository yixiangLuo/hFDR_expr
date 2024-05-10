# Numerical experiments for hFDR

## Description and structure

This repo hosts the code used to generate figures in the cKnockoff paper.

`R/`: host R code for the numerical experiments

`data/`: store the numerical results

`figs/`: `R/plot.R` will read the numerical results in `data/` and save the figures in `figs/`

`jobs/`: Slurm job files to be submitted to cluster for large scale simulation

`log/`: store the Slurm logs.

### Dependency

Required packages:
- hFDR: `devtools::install_github("yixiangluo/hFDR")`
- glmnet
- glasso
- CVXR
- here
- abind
- foreach
- doParallel

## Usage

### Simulation studies

#### local experiments

1. Create a `my_expr_setting.R` file in `R/settings` (see existed examples), in which specify the experiment settings.
2. Run `Rscript R/do_expr.R my_expr_setting` then `draw_figs.R` with corresponding settings.
3. The real time progress are shown in `data/temp/progress-my_expr_setting.txt`, the results are stored in `data/my_expr_setting.RData`, and the produced figure is in `figs/`.

#### cluster experiments

1. Echo the names of the experiments (e.g. `my_expr_setting`) in `jobs/expr_names.sh`
2. Specify the Slurm arguments in `jobs/expr_template.sh` according your cluster settings. e.g. change the notification email at line 11. But please do NOT add/delete lines or change the order of lines, since this is a template file the content will be read and modified by `R/cluster/set_expr.R` according to which line the content is in.
3. Run `jobs/set_exprs.sh` (e.g. `sbatch jobs/set_exprs.sh`). This will prepare the shared data and structure used by each task in cluster computing.
4. Submit the jobs to the cluster by `sh jobs/do_exprs.sh`.
5. Check the progress at `sh R/cluster/show_progress.R`.
6. After all the experiments are done, run `sh jobs/post_process.sh`. This will aggregate the numerical results from all tasks and draw figures.
7. If there is an error, check the log in `log/`.

### Illustrative figures in the paper

To replicate, please run `R/illustration/illustration.R`.

### HIV drug resistance

To replicate, please run the files in `R/HIV/HIV_preprocess.R` and then `R/HIV/HIV.R`.

### Protein signaling network

To replicate, please run `R/protein/protein.R`.
