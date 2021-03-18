# Pipeline containers

To have the required containers available, they must be built locally. For distribution license reasons this repository does not contain the containers itself but the recipes to easily build locally (requires [singularity](https://sylabs.io/docs/)). Examples can be found in [build_local.sh](build_local.sh).
For the containers found on DockerHub, it is simply:
```
singularity build upstream/openswath:0.1.2.simg docker://openswath/openswath:0.1.2
singularity build upstream/pyprophet_2.1.3.simg docker://pyprophet/pyprophet:2.1.3
```
For a specific container that is not on DockerHub, change into the directories for the chosen container that contains a `*.sdef` file. From there, it is for example:
```
singularity build downstream/pridereswath-downstream.simg python-r:bioc-py3-tinytex-qc_plus_msstats.sdef
```
__Note:__ After containers are built, the workflow [parameter file](../inputs/configs/pridereswath_v1.1.yaml) needs to be adjusted to reflect the location of the respective containers.


