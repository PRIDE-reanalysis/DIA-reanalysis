#!/usr/bin/env bash

# upstream
singularity build upstream/openswath:0.1.2.simg docker://openswath/openswath:0.1.2
singularity build upstream/pyprophet_2.1.3.simg docker://pyprophet/pyprophet:2.1.3
singularity build upstream/sciex_converter_prideswath.simg container/upstream/wiffconverter/sciex_converter_prideswath.sdef
singularity build upstream/pridereswath.simg container/upstream/pridereswath/pridereswath-upstream.sdef 
singularity build upstream/yamato:v1.0.4_pridereswath.simg container/upstream/yamato:v1.0.4_pridereswath.sdef

# downstream
# from docker builds:
#sudo docker build -t python-r:bioc-py3-tinytex-qc_plus_msstats .
#sudo singularity build python-r:bioc-py3-tinytex-qc_plus_msstats.simg docker-daemon://python-r:bioc-py3-tinytex-qc_plus_msstats
# direct builds:
singularity build downstream/pridereswath-downstream.simg downstream/python-r:bioc-py3-tinytex-qc_plus_msstats.sdef

# postprocess
# see downstream
