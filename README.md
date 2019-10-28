A Nextflow pipeline for scmap package. The workflow accepts expression values for query and reference datasets as 10X directories. It then runs the necessary pre-processing, feature selection and indexing steps. Depending on parameter specified in `nextflow.config`, the workflow produces either cell-to-cluster or cell-to-cell projections. Schematic diagram below shows the workflow structure. 
![](https://github.com/ebi-gene-expression-group/scmap-workflow/blob/develop/scmap_nextflow.png)   