# Overview
The general workflow consists of these steps:
(0. Curate experimental design, wiff/scan files, iRT and assay librarary targets)
1. Converting the RAW data into mzML file
2. Extract window setting, QC and optimise assay lib
   - check window consistency
   - create SwaMe mzQC
   - generate decoys and optimise assay lib 
3. OpenSWATH
4. Pyprophet 
   - Merge `.osw` files
   - Training a model (with subsampling and backpropagation)
   - Apply model for FDR control 
5. TRIC alignment
(6. Statistical analysis)

## The general workflow is split into an automated part and an supervised part
The general workflow is split into an automated and an supervised part. Steps 1 to 5 are all included in the nextflow workflow. Steps 0 and 6 need to be executed manually. (To find out how [see here for input preparation](lsf.md) and [here for MSStats input](http://openswath.org/en/latest/docs/swath2stats.html).)

For obvious reasons, there needs to be some 'information' collected and more importantly the data input organised accordingly prior to any automated workflow excution. However, the experimental design does not influence how the SWATH runs are analysed individually. Thus, experimental design and MSstats script customisation can be created during nextflow excution. 

### Automated workflow
To meet the technical requirements of nextflow, the general workflow steps are arranged like this.

![nextflow](flowchart_edited.svg "The visualised nextflow")
> generated with `nextflow_20.01.0 run ../workflows/openswath-nextflow_0.6.nf -with-dag flowchart.dot` and corrections in the dot file to `dot -Tpng flow.dot > dot.png`

#### nextflow parameters
There are 5 parameters for running the nextflow, for which execution documentation can also be had with `nextflow_20.01.0 run openswath-nextflow-no_scaleup.nf --help`:
  - `--wiff_ftp` The PRIDE ftp raw file URI
  - `--scan_ftp` The PRIDE ftp raw file URI
  - `--lib` The assay library
  - `--irtkit` The traml with iRT peptides
  - `--out_dir` The (existing) directory to deposit the TRIC results (and the merged OSW file)

The `--wiff_ftp` and `--scan_ftp` parameters should point to the 'directory' that holds the sciex raw file pairs (`.wiff`/`.wiff.scan` pairs). 
These are handled as URIs which means either a URL to a FTP/HTTP directory or a linux filesystem path. 
In the latter case, the directories may contain symlinks, so no copy of the files is necessary. 
It is however necessary, that the filepairs are named the *same*, 
except their file extension which should be either `.wiff` or `.wiff.scan` 
(*not* `.wiff`/`.scan` - this is a hardcoded requirement of the sciex converter).  
It is okay if both parameters have the same value, but not if one is empty or omitted.

The `--lib` parameter should always point with an absolute path to the `.tsv` file that fits the experiment files. 
It should contain all targets and decoy. Preparation steps are documented in [the queries documentation](documentation/assaylib.md)

The `--irtkit` parameter specifies which iRT peptides were used in the experiment's runs and therfore should be considered for guiding the analysis. 
The most used 'kit' is the BiognoSys iRT kit which is added to the sample just before the MS run.
The file type should be TraML, as described [here](http://biorxiv.org/content/early/2016/03/19/044552).

The last parmeter is `--out_dir` which should point to the already existing folder into which the TRIC results (and the merged OSW file) will be deposited.

TODO: Evaluate if specifying the analysis SWATH windows is still necessary.


#### Implicit behaviour
Each assay library is genereated from scratch. It is optimised for the current window setting extended with decoy addition. Details [are here](documentation/assaylib.md).

