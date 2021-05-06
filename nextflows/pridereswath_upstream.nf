#!/usr/bin/env nextflow

/*
===============================================================================
            PRIDE re-SWATH upstream workflow v1.0.0
===============================================================================
 @authors
 Mathias Walzer <walzer@ebi.ac.uk>
-------------------------------------------------------------------------------
Pipeline overview:
 - 0:   Converting the RAW data into mzML file
 - 1:   Extract window setting and create/optimise library-assay
 - 2:   OpenSWATH analysis
 - 3:   Pyprophet model merging, training, scoring, and export
 - 4:   TRIC alignment
-------------------------------------------------------------------------------
Some notes on the technical aspects of the workflow:
* The workflow is designed for use with HPC and containers.
* To start, you should have this file ready together with the `params-file` and
    a `nextflow.config` file customised to your `cluster`. The `params-file`
    needs customisation regarding your choice of container use. There, you can
    also adjust workflow parameters to the individual process steps.
* The workflow depends on the correct container configuration for each process
    step, to be defined in the `params-file`.
* Due to the generous sizes of mzMLs converted from wiff/scans, we need to use
    disk caching in OpenSwathWorkflow. Even then, because of the size, we might
    not be able to rely on `/scratch` (not present on all nodes or full 
    already) - therefore the cache directory is fixed to the output directory, 
    but is cleaned after each successful process finish. Be sure to designate 
    your output directory to a location with adequate I/O specs (e.g. not a 
    backup device)
* In case the workflow execution fails or nextflow does not clean up after
    itself even if the parameters are set in the nextflow.config (NB `remove`
    is just for the container) you need to manually remove the tmp_cache and
    work folders in the output directory.
* The success of using container URLs in the configuration instead of a
    singularity image path may depend on the container integration of the
    cluster configuration. You should probably add the cacheDir to your
    nextflow.config
* The output of the workflow is the created assay libray (pqp), the individual
    runs' window setting (txt), mzQC files from SwaMe, OpenSwathWorkflow
    result files, and final TRIC results (tsv).
*/

def helpMessage() {
    log.info"""
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run pridereswath.nf -c nextflow.config -profile local -params-file pridereswath.yaml <arguments>

    Workflow arguments:

    --wiff The dataset location (Can be symlinks,
                only wiff/scan pairs are considered)
    --lib The library file (in osw compatible tsv format, without decoys)
    --no_add_decoy If set, workflow will skip decoy addition (make sure it
                contains decoys; the file to `--lib parameter` can now also be
                in pqp or TraML format)
    --irtkit The TraML file containing the iRT peptides chosen (recommended)
    --out_dir results folder
    
    Nextflow general arguments (see nextflow documentation):
    -params-file A YAML configuration file for workflow parameters
                (note the single dash for the parameter name)
    -c A nextflow config file
    -profile compute profile as defined in the nextflow config file


    An example call on LSF would look like this:
    ```
bsub -J test -M 4096 -R "rusage[mem=4096]" /nfs/nextflow_20.01.0 run /nfs/workflows/pridereswath.nf -params-file /nfs/test/input/pridereswath.yaml -c /nfs/workflows/nextflow.config -profile cluster --irtkit /nfs/test/ion_libraries/irt-kit-reference-sheet.TraML --lib /nfs/test/ion_libraries/phl004_canonical_sall_osw.tsv --wiff /nfs/test/input/ --out_dir /nfs/test/output/
    ```
    """.stripIndent()
}


/*
 * Config Variables
 */
// show help message if requested
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// defaults where conditionaly needed
params.wiff = false
params.no_add_decoy = false
params.sample_size = 0
params.direct_osw = false
params.direct_pqp = false

// prepare input (into value & channels)
if ( (params.direct_osw && !params.direct_pqp) || (!params.direct_osw && params.direct_pqp) ){
    println("direct mode only with osw AND pqp input possible!")
    exit 0
}

if (params.direct_osw && params.direct_pqp){
    println("using direct osw feed!")
    direct_osw_1 = Channel.fromPath( "${params.direct_osw}/*.osw" )
    direct_osw_2 = Channel.fromPath( "${params.direct_osw}/*.osw" )
    //osw_channel_1 = direct_osw_1
    //osw_channel_2 = direct_osw_2
    println("using direct pqp feed!")
    //pqp_channel = direct_pqp
    direct_pqp = Channel.fromPath( "$params.direct_pqp" )
    wiffFiles = Channel.empty()
} else {
    wiffFiles = Channel.fromFilePairs( "${params.wiff}/*{.wiff,.wiff.scan}" ) 
    direct_osw_1 = Channel.empty()
    direct_osw_2 = Channel.empty()
    direct_pqp = Channel.empty()
}

// calculate subsample_ratio
params.sample_size.println()
if (params.sample_size.intValue() <= 0) {
    println('re-eval sample size from channel size')
    if ( params.direct_osw ){
        tmpchannel = Channel.fromPath( "${params.direct_osw}/*.osw" )
        sample_size = tmpchannel.count().floatValue().value
    }else{
        tmpchannel = Channel.fromFilePairs( "${params.wiff}/*{.wiff,.wiff.scan}" )
        sample_size = tmpchannel.count().floatValue().value
    }
    //sample_size = tmpchannel.count().floatValue().value
    println "overriding sample size $sample_size"
}
else {
    sample_size = params.sample_size.floatValue()
}
println "input sample size: $sample_size"
subsample_ratio = 1.0 / sample_size.floatValue()
println "subsample ratio: $subsample_ratio"
// println "by .value ${subsample_ratio.value}"


/*
 * Process definitions
 */
process wiffCONVERT {
    container "${params.sciex_converter.container}"
    memory { params.sciex_converter.memory.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries { params.sciex_converter.maxRetries }

    input:
    tuple val(wiff_base), file(filepaths) from wiffFiles

    output:
    file '*.mzML' into mzmlwindow

    when:
    ( params.wiff )

    // --index does not work with converter version 0.7, versions above do not work properly
    """
    /convert_wiff_ignore_NullReferenceException.sh ${wiff_base}.wiff ${wiff_base}.mzML
    """
}


process mzmlWINDOW {
    container "${params.helper_scripts.container}"
    memory { params.helper_scripts.memory.GB * task.attempt }
    errorStrategy 'retry'
    publishDir "${params.out_dir}/${mzML_file.baseName}/" , mode: 'copy', pattern: "*.txt"

    input:
    file mzML_file from mzmlwindow

    output:
    tuple file(mzML_file), file('*_windows.txt') into mzml_channel
    file('*_windows.txt') into checkwindows
    file('*_windows.txt') into window4pqp_channel
    file(mzML_file) into checkqc

    when:
    ( params.wiff )

    script:
    """
    python3 /extract_windows.py --in ${mzML_file}
    """
}


process windowCHECK {
    container "${params.helper_scripts.container}"
    memory { params.helper_scripts.memory.GB * task.attempt }
    errorStrategy 'retry'

    input:
    file(p) from checkwindows.collect()

    when:
    ( params.wiff )

    script:
    """
    python3 /compare_windows.py --in ${p}
    """
}


process mzmlQC {
    container "${params.yamato.container}"
    memory { params.yamato.memory.GB * task.attempt }
    errorStrategy { task.exitStatus in 0..1 ? 'ignore' : 'retry' }
    maxRetries { params.yamato.maxRetries }
    cpus { params.openswath.threads }
    publishDir "${params.out_dir}/${mzML_file.baseName}/" , mode: 'copy', pattern: "*.mzQC"

    input:
    file(mzML_file) from checkqc

    output:
    file('*.mzQC') into mzQC_channel

    when:
    ( params.wiff )

    script:
    if ( params.irtkit )
        """
        mkdir -p ${params.out_dir}/tmp_cache/${mzML_file.baseName}_SwaMe
        /swame_ignore_throws.sh ${mzML_file} ${mzML_file.baseName}.mzQC $params.irtkit ${params.out_dir}/tmp_cache/${mzML_file.baseName}_SwaMe
        rm -rf ${params.out_dir}/tmp_cache/${mzML_file.baseName}_SwaMe
        """
    else
        """
        mkdir -p ${params.out_dir}/tmp_cache/${mzML_file.baseName}_SwaMe
        /swame_ignore_throws.sh ${mzML_file} ${mzML_file.baseName}.mzQC $params.irtkit ${params.out_dir}/tmp_cache/${mzML_file.baseName}_SwaMe
        rm -rf ${params.out_dir}/tmp_cache/${mzML_file.baseName}_SwaMe
        """
}


process libOSW {
    container "${params.openswath.container}"
    memory { params.openswath.memory.GB * task.attempt }
    errorStrategy 'retry'
    cpus { params.openswath.threads }
    publishDir "${params.out_dir}/" , mode: 'copy', pattern: "*.pqp"

    input:
    file(custom_windows_file) from window4pqp_channel.first()

    output:
    file('*.pqp') into pqp_channel

    when:
    ( params.wiff )

    script:
    if ( params.no_add_decoy )
        """
        TargetedFileConverter -in $params.lib -out assay_plus_decoys.pqp -threads ${task.cpus}
        """
    else
        """
        OpenSwathAssayGenerator \
        -threads ${task.cpus} \
        -min_transitions 3 \
        -max_transitions 6 \
        -allowed_fragment_types b,y \
        -allowed_fragment_charges 1,2,3,4 \
        -enable_detection_specific_losses \
        -enable_detection_unspecific_losses \
        -precursor_mz_threshold 0.025 \
        -precursor_lower_mz_limit 400 \
        -precursor_upper_mz_limit 1200 \
        -product_mz_threshold 0.025 \
        -product_lower_mz_limit 350 \
        -product_upper_mz_limit 2000 \
        -in $params.lib \
        -out assay.TraML \
        -swath_windows_file ${custom_windows_file}

        OpenSwathDecoyGenerator \
        -threads ${task.cpus} \
        -allowed_fragment_types b,y \
        -allowed_fragment_charges 1,2,3,4 \
        -enable_detection_specific_losses \
        -enable_detection_unspecific_losses \
        -product_mz_threshold 0.025 \
        -in assay.TraML \
        -out assay_decoys.TraML

        TargetedFileConverter -in assay_decoys.TraML -out assay_plus_decoys.pqp -threads ${task.cpus}
        """
}


process mzmlOSW {
    container "${params.openswath.container}"
    memory { params.openswath.memory.GB * 2 * task.attempt }
    errorStrategy 'retry'
    cpus { params.openswath.threads }
    publishDir "${params.out_dir}/${mzML_file.baseName}/" , mode: 'copy', pattern: "*.osw", overwrite: "true"

    input:
    tuple file(mzML_file), file(custom_windows_file) from mzml_channel
    file(pqp_file) from pqp_channel
    
    output:
    file('*.osw') into osw_channel_1
    file('*.osw') into osw_channel_2

    when:
    ( params.wiff )

    script:
    if ( params.irtkit )
        """
        mkdir -p ${params.out_dir}/tmp_cache/${mzML_file.baseName}
        OpenSwathWorkflow \
            -threads ${task.cpus} \
            -in ${mzML_file} \
            -tr ${pqp_file} \
            -tr_irt $params.irtkit \
            -swath_windows_file ${custom_windows_file} \
            -out_osw ${mzML_file.baseName}.osw \
            -min_upper_edge_dist 0 \
            -batchSize 1000 \
            -force \
            -use_ms1_traces \
            -Scoring:Scores:use_ms1_mi \
            -Scoring:Scores:use_mi_score \
            -readOptions cacheWorkingInMemory \
            -tempDirectory ${params.out_dir}/tmp_cache/${mzML_file.baseName}/ \
            -mz_extraction_window $params.openswath.mz_extraction_window \
            -ppm \
            -mz_correction_function unweighted_regression \
            -irt_mz_extraction_window $params.openswath.irt_mz_extraction_window \
            -ppm_irtwindow \
            -rt_extraction_window $params.openswath.rt_extraction_window \
            -Scoring:stop_report_after_feature 5 \
            -Scoring:TransitionGroupPicker:compute_peak_quality true \
            -Scoring:TransitionGroupPicker:compute_peak_shape_metrics \
            -Scoring:TransitionGroupPicker:background_subtraction ${params.openswath.background_subtraction} \
            -Scoring:TransitionGroupPicker:min_peak_width ${params.openswath.min_peak_width}
        rm -rf ${params.out_dir}/tmp_cache/${mzML_file.baseName}
        """
    else
        """
        mkdir -p ${params.out_dir}/tmp_cache/${mzML_file.baseName}
        OpenSwathWorkflow \
            -threads ${task.cpus} \
            -in ${mzML_file} \
            -tr ${pqp_file} \
            -swath_windows_file ${custom_windows_file} \
            -out_osw ${mzML_file.baseName}.osw \
            -min_upper_edge_dist 0 \
            -batchSize 1000 \
            -force \
            -use_ms1_traces \
            -Scoring:Scores:use_ms1_mi \
            -Scoring:Scores:use_mi_score \
            -readOptions cacheWorkingInMemory \
            -tempDirectory ${params.out_dir}/tmp_cache/${mzML_file.baseName}/ \
            -mz_extraction_window $params.openswath.mz_extraction_window \
            -ppm \
            -mz_correction_function unweighted_regression \
            -irt_mz_extraction_window $params.openswath.irt_mz_extraction_window \
            -ppm_irtwindow \
            -rt_extraction_window $params.openswath.rt_extraction_window \
            -Scoring:stop_report_after_feature 5 \
            -Scoring:TransitionGroupPicker:compute_peak_quality true \
            -Scoring:TransitionGroupPicker:compute_peak_shape_metrics \
            -Scoring:TransitionGroupPicker:background_subtraction ${params.openswath.background_subtraction} \
            -Scoring:TransitionGroupPicker:min_peak_width ${params.openswath.min_peak_width}
        rm -rf ${params.out_dir}/tmp_cache/${mzML_file.baseName}
        """
        // 2.5 only: Unknown option(s) '[x_ -ms1_isotopes, r_ -irt_mz_extraction_window_unit, r_ -mz_extraction_window_ms1_unit, x_ -mz_extraction_window_ms1, x_ -mz_extraction_window_unit]'  // x: deleted r: replacement
}


process subsamplePYPROPHET {
    container "${params.pyprophet.container}"
    memory { params.pyprophet.memory.GB * task.attempt }
    errorStrategy 'retry'

    input:
    file(osw_file) from osw_channel_1.concat( direct_osw_1 )

    output:
    file('*.osws') into subsample_channel

    script:
    if ( params.pyprophet.test )
        """
        pyprophet subsample --in=${osw_file} --out=${osw_file.baseName}.osws --subsample_ratio=${subsample_ratio.value} --test
        """
    else
        """
        pyprophet subsample --in=${osw_file} --out=${osw_file.baseName}.osws --subsample_ratio=${subsample_ratio.value}
        """
}


process modelPYPROPHET {
    container "${params.pyprophet.container}"
    memory { params.pyprophet.memory.GB * task.attempt }
    errorStrategy 'retry'
    cpus 1

    input:
    file(pqp_file) from pqp_channel.concat( direct_pqp )    
    file('*.osws') from subsample_channel.collect()

    output:
    file 'model.oswm' into subsamplemodel_channel_1
    file 'model.oswm' into subsamplemodel_channel_2

    script:
    if ( params.pyprophet.test )
        """
        pyprophet merge --template=${pqp_file} --out=model.oswm *.osws
        pyprophet score --in=model.oswm --level=ms1ms2 --tric_chromprob --test
        """
    else
        """
        pyprophet merge --template=${pqp_file} --out=model.oswm *.osws
        pyprophet score --in=model.oswm --level=ms1ms2 --tric_chromprob
        """
}


process applyPYPROPHET {
    container "${params.pyprophet.container}"
    memory { params.pyprophet.memory.GB * task.attempt }
    errorStrategy 'retry'

    input:
    file(model_file) from subsamplemodel_channel_1.first()
    file(osw_file) from osw_channel_2.concat( direct_osw_2 )
    // tuple file(model_file), file(osw_file) from subsamplemodel_channel_1.combine( direct_osw_2 )
    // file(osw_file) from osw_channel_2

    output:
    file('*.oswr') into rescored_channel
    file(osw_file) into osw_channel_3

    script:
    if ( params.pyprophet.test )
        """
        pyprophet score --in=${osw_file} --apply_weights=${model_file} --level=ms1ms2 --tric_chromprob --test
        pyprophet reduce --in=${osw_file} --out=${osw_file.baseName}.oswr
        """
    else
        """
        pyprophet score --in=${osw_file} --apply_weights=${model_file} --level=ms1ms2 --tric_chromprob
        pyprophet reduce --in=${osw_file} --out=${osw_file.baseName}.oswr
        """

}


process mergePYPROPHET {
    container "${params.pyprophet.container}"
    memory { params.pyprophet.memory.GB * task.attempt }
    errorStrategy 'retry'

    input:
    file(model_file) from subsamplemodel_channel_2.first()
    file(osw_file) from rescored_channel.collect()

    output:
    file('model_global.osw') into globalmodel_channel

    """
    pyprophet merge --template=${model_file} --out=model_global.osw *.oswr
    pyprophet peptide --context=global --in=model_global.osw
    pyprophet protein --context=global --in=model_global.osw
    """

}


process backpropPYPROPHET {
    container "${params.pyprophet.container}"
    memory { params.pyprophet.memory.GB * task.attempt }
    errorStrategy 'retry'
    publishDir "${params.out_dir}/pyprophet_export/" , mode: 'copy', pattern: "*.tsv", overwrite: "true"

    input:
    // file(osw_file) from osw_channel_3
    // file(model_global) from globalmodel
    tuple file(osw_file), file(model_global) from osw_channel_3.combine( globalmodel_channel )

    output:
    file('*_export.tsv') into pyprophecies_channel

    """
    pyprophet backpropagate --in=${osw_file} --apply_scores=${model_global}
    pyprophet export --in=${osw_file} --out=${osw_file.baseName}_export.tsv \
        --format=legacy_merged \
        --max_rs_peakgroup_qvalue $params.pyprophet.max_rs_peakgroup_qvalue \
        --max_global_peptide_qvalue $params.pyprophet.max_global_peptide_qvalue \
        --max_global_protein_qvalue $params.pyprophet.max_global_peptide_qvalue
    """
}


process alignTRIC {
    container "${params.tric.container}"
    memory { params.tric.memory.GB * task.attempt }
    errorStrategy 'retry'
    cpus { params.tric.threads }
    publishDir params.out_dir, mode: 'move'

    input:
    file(merged_exports) from pyprophecies_channel.collect()

    output:
    file 'tric*.tsv' into result

    """
    feature_alignment.py \
        --in *.tsv \
        --out tric.tsv \
        --out_matrix tric_matrix.tsv \
        --method LocalMST \
        --realign_method lowess \
        --max_rt_diff 60 \
        --mst:useRTCorrection True \
        --mst:Stdev_multiplier 3.0 \
        --target_fdr -1 \
        --fdr_cutoff $params.tric.fdr_cutoff \
        --max_fdr_quality $params.tric.max_fdr_quality \
        --alignment_score $params.tric.alignment_score
    """
}
