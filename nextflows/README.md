# PRIDEreSWATH

This folder contains the workflow design files for the nextflow workflows.
The respective main workflow files are the `*.nf` files, however you may need additional input and configuration files which might also need adjustment to your infrastructure, found in the `inputs` folder, for more details see [below](## Documentation). Workflow execution will also depend on the availability of the containers, read on in the documentation and [below](## Documentation). 


## Documentation
The documentation helpful for technical aspects of the workflow execution are in the `docs` folder.

## Container
Singularity container definition files referenced in the workflow configuration file can be found in the `container` folder. Brief build instructions can also be found there, for more see the documentation.

## Quickstart
```
./nextflow run pridereswath.nf -params-file pridereswath.yaml -c nf.config --irtkit /fs/irt.TraML --lib /fs/library.tsv --wiff /fs/input/ --out_dir /fs/output/
```
