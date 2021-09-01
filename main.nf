#! usr/bin/env nextflow

nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
reps = 10000

download_vcf = null
params.R_libpath = ""