#! usr/bin/env nextflow

nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )
reps = 10000

vcf_file=chanel.fromPath("/projects/b1059/projects/Nicolas/c.briggsae/variant_calling/WI-20210803/variation/WI.20210803.hard-filter.vcf.gz")
vcf_index=channel.fromPath("/projects/b1059/projects/Nicolas/c.briggsae/variant_calling/WI-20210803/variation/WI.20210803.hard-filter.vcf.gz.tbi")
params.traitfile("heritability calc data path")


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