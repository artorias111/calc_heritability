#! usr/bin/env nextflow

// if( !nextflow.version.matches('20.0+') ) {
//     println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
//     println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
//     exit 1
// } this is a test

nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
reps = 10000


/*
~ ~ ~ > * Parameters
*/

download_vcf = null
params.R_libpath = ""

// VCF param
if(params.debug) {
    println """
        *** Using debug mode ***
    """
    // debug for now with small vcf
    params.vcf = "330_TEST.vcf.gz"

    vcf_file = Channel.fromPath("${binDir}/test_data/330_TEST.vcf.gz")
    vcf_index = Channel.fromPath("${binDir}/test_data/330_TEST.vcf.gz.tbi")
    params.traitfile = "${binDir}/test_data/ExampleTraitData.csv"

    // lower number of reps for debug
    reps = 100
        
} else if(params.gcp) { 
    // use the data directly from google on gcp
    vcf_file = "gs://elegansvariation.org/releases/20210121/variation/WI.20210121.hard-filter.isotype.vcf.gz"
    vcf_index = "gs://elegansvariation.org/releases/20210121/variation/WI.20210121.hard-filter.isotype.vcf.gz.tbi"

} else if(!params.vcf) {
    // if there is no VCF date provided, pull the latest vcf from cendr.
    params.vcf = "20210121"
    download_vcf = true
    
} else {
    // use the vcf data from QUEST when a cendr date is provided

    // Check that params.vcf is valid
    if("${params.vcf}" == "20210121" || "${params.vcf}" == "20200815" || "${params.vcf}" == "20180527" || "${params.vcf}" == "20170531") {
        vcf_file = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.hard-filter.isotype.vcf.gz")
        vcf_index = Channel.fromPath("/projects/b1059/data/c_elegans/WI/variation/${params.vcf}/vcf/WI.${params.vcf}.hard-filter.isotype.vcf.gz.tbi")

    } else {
        println "Error! Cannot find ${params.vcf}.vcf. Please provide a valid CeNDR release date (20210121, 20200815, 20180527, or 20170531)."
        exit 1
    }
}



/*
~ ~ ~ > * WORKFLOW
*/

workflow {

	// if no VCF is provided, download the latest version from CeNDR
	if(download_vcf) {
	    pull_vcf()

	    vcf_file = pull_vcf.out.hard_vcf
	    vcf_index = pull_vcf.out.hard_vcf_index
	}

	// Fix strain names
    Channel.fromPath("${params.traitfile}")
        .combine(Channel.fromPath("${binDir}/bin/strain_isotype_lookup.tsv"))
        .combine(Channel.fromPath("${binDir}/bin/Fix_Isotype_names_bulk_h2.R")) | fix_strain_names_bulk        

    traits_to_map = fix_strain_names_bulk.out.fixed_strain_phenotypes
            .flatten()
            .map { file -> tuple(file.baseName.replaceAll(/pr_/,""), file) }

    // Genotype matrix
    pheno_strains = fix_strain_names_bulk.out.phenotyped_strains_to_analyze

    vcf_file.spread(vcf_index)
            .combine(pheno_strains) | vcf_to_geno_matrix

    // calclate heritability and generate report output
    traits_to_map 
    	.combine(vcf_to_geno_matrix.out)
        .combine(Channel.from(reps))
        .combine(Channel.fromPath("${binDir}/bin/20210716_H2_script.R")) | heritability

    //generate html report
    traits_to_map
    	.join(heritability.out)
    	.combine(fix_strain_names_bulk.out.strain_issues)
        .combine(Channel.fromPath("${binDir}/bin/20210716_hert_report.Rmd")) | html_report

}

/*
==============================================
~ > *                                    * < ~
~ ~ > *                                * < ~ ~
~ ~ ~ > *  DOWNLOAD VCF FROM CENDR   * < ~ ~ ~
~ ~ > *                                * < ~ ~
~ > *                                    * < ~
==============================================
*/

process pull_vcf {

    tag {"PULLING VCF FROM CeNDR"}

    output:
        path "*hard-filter.isotype.vcf.gz", emit: hard_vcf 
        path "*hard-filter.isotype.vcf.gz.tbi", emit: hard_vcf_index 
        path "*impute.isotype.vcf.gz", emit: impute_vcf 
        path "*impute.isotype.vcf.gz.tbi", emit: impute_vcf_index 
        path "*.strain-annotation*.tsv", emit: ann_vcf

    """
        wget https://storage.googleapis.com/elegansvariation.org/releases/${params.vcf}/variation/WI.${params.vcf}.small.hard-filter.isotype.vcf.gz
        tabix -p vcf WI.${params.vcf}.small.hard-filter.isotype.vcf.gz

        wget https://storage.googleapis.com/elegansvariation.org/releases/${params.vcf}/variation/WI.${params.vcf}.impute.isotype.vcf.gz
        tabix -p vcf WI.${params.vcf}.impute.isotype.vcf.gz

        wget https://storage.googleapis.com/elegansvariation.org/releases/${params.vcf}/variation/WI.${params.vcf}.strain-annotation.bcsq.tsv

    """
}

/*
==============================================================
~ > *                                                    * < ~
~ ~ > *                                                * < ~ ~
~ ~ ~ > *  FIX STRAIN NAMES TO MATCH THOSE ON CENDR  * < ~ ~ ~
~ ~ > *                                                * < ~ ~
~ > *                                                    * < ~
==============================================================
*/

/*
THIS WILL NEED TO BE UPDATED TO HANDLE OTHER SPECIES
*/


process fix_strain_names_bulk {

    tag {"BULK TRAIT"}

    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "*pr_*.tsv"
    publishDir "${params.out}/Phenotypes", mode: 'copy', pattern: "strain_issues.txt"

    input:
        tuple file(phenotypes), file(isotype_lookup), file(fix_isotype_script)

    output:
        path "pr_*.tsv", emit: fixed_strain_phenotypes 
        path "Phenotyped_Strains.txt", emit: phenotyped_strains_to_analyze 
        path "strain_issues.txt", emit: strain_issues

    """
        # add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${fix_isotype_script} > Fix_Isotype_names_bulk.R 
        Rscript --vanilla Fix_Isotype_names_bulk.R ${phenotypes} fix $isotype_lookup
    """

}

/*
===========================================================
~ > *                                                 * < ~
~ ~ > *                                             * < ~ ~
~ ~ ~ > *  CONVERT THE VCF TO A GENOTYPE MATRIX   * < ~ ~ ~
~ ~ > *                                             * < ~ ~
~ > *                                                 * < ~
===========================================================
*/

process vcf_to_geno_matrix {

    //publishDir "${params.out}/Genotype_Matrix", mode: 'copy'

    cpus 1

    input:
        tuple file(vcf), file(index), file(strains)

    output:
        file("Genotype_Matrix.tsv") 

    """
        bcftools view -S ${strains} ${vcf} |\\
        bcftools filter -i N_MISSING=0 -Oz -o Phenotyped_Strain_VCF.vcf.gz
        tabix -p vcf Phenotyped_Strain_VCF.vcf.gz
        plink --vcf Phenotyped_Strain_VCF.vcf.gz \\
            --snps-only \\
            --biallelic-only \\
            --maf 0.05 \\
            --set-missing-var-ids @:# \\
            --indep-pairwise 50 10 0.8 \\
            --geno \\
            --allow-extra-chr
        awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
        sort -k1,1d -k2,2n > markers.txt
        bcftools query -l Phenotyped_Strain_VCF.vcf.gz |\\
        sort > sorted_samples.txt 
        bcftools view -v snps \\
        -S sorted_samples.txt \\
        -R markers.txt \\
        Phenotyped_Strain_VCF.vcf.gz |\\
        bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' |\\
        sed 's/[[# 0-9]*]//g' |\\
        sed 's/:GT//g' |\\
        sed 's/0|0/-1/g' |\\
        sed 's/1|1/1/g' |\\
        sed 's/0|1/NA/g' |\\
        sed 's/1|0/NA/g' |\\
        sed 's/.|./NA/g'  |\\
        sed 's/0\\/0/-1/g' |\\
        sed 's/1\\/1/1/g'  |\\
        sed 's/0\\/1/NA/g' |\\
        sed 's/1\\/0/NA/g' |\\
        sed 's/.\\/./NA/g' > Genotype_Matrix.tsv
    """

}


/*
===========================================================
~ > *                                                 * < ~
~ ~ > *                                             * < ~ ~
~ ~ ~ > * CALCULATE BROAD AND NARROW HERITABILITY * < ~ ~ ~
~ ~ > *                                             * < ~ ~
~ > *                                                 * < ~
===========================================================
*/

process heritability {

	publishDir "${params.out}/", mode: 'copy'

	input:
		tuple val(TRAIT), file(phenotype), file(geno_matrix), val(reps), file(h2_script)


	output:
		tuple val(TRAIT), file("heritability_result.tsv")


	"""
		# add R_libpath to .libPaths() into the R script, create a copy into the NF working directory 
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${h2_script} > H2_script.R 
        Rscript --vanilla H2_script.R ${phenotype} ${geno_matrix} ${reps}

	"""


}


/*
===========================================================
~ > *                                                 * < ~
~ ~ > *                                             * < ~ ~
~ ~ ~ > *           GENERATE HTML REPORT          * < ~ ~ ~
~ ~ > *                                             * < ~ ~
~ > *                                                 * < ~
===========================================================
*/

process html_report {

	publishDir "${params.out}/", mode: 'copy'

	input:
		tuple val(TRAIT), file(phenotype), file(hert), file(strain_issues), file(html_report)

	output:
		tuple file("*.html"), file("h2_plot.png")


	"""
		cat "${html_report}" | \\
		sed 's+Phenotypes/pr_TRAITNAME.tsv+${phenotype}+g' | \\
		sed "s+TRAITNAME+${TRAIT}+g" | \\
		sed 's+heritability_result.tsv+${hert}+g' | \\
		sed 's+Phenotypes/strain_issues.txt+${strain_issues}+g' > hert_report_${TRAIT}.Rmd 
	    echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" > .Rprofile
	    Rscript -e "rmarkdown::render('hert_report_${TRAIT}.Rmd')"

	"""


}






