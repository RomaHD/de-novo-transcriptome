export NXF_HOME=/omics/groups/OE0436/data/kurilov/nikolas/oxford_nanopore/nf_cache

export NXF_CONDA_CACHEDIR=${NXF_HOME}/conda

export NXF_SINGULARITY_CACHEDIR=${NXF_HOME}/apptainer

export NXF_DISABLE_CHECK_LATEST=true

export SINGULARITY_DISABLE_CACHE=True

module load nextflow/24.10.2
nextflow run epi2me-labs/wf-transcriptomes -profile singularity --fastq /omics/odcf/project/OE0420/setdb1_human/nonOTP/ont/raw/i0043665/data/1616354/fastq_pass --out_dir /omics/groups/OE0436/data/kurilov/nikolas/oxford_nanopore --cdna_kit SQK-PCS114 --ref_annotation /omics/groups/OE0436/data/kurilov/nikolas/de_novo_colo800/gffcompare_output/gffCompare.annotated.sorted_ver2.gtf --ref_genome /omics/groups/OE0436/data/kurilov/nikolas/de_novo_rna2/reference_data/hg19.fa --sample COLO800
nextflow run epi2me-labs/wf-transcriptomes -profile singularity --fastq /omics/odcf/project/OE0420/setdb1_human/nonOTP/ont/raw/i0043665/data/1616356/fastq_pass --out_dir /omics/groups/OE0436/data/kurilov/nikolas/oxford_nanopore --cdna_kit SQK-PCS114 --ref_annotation /omics/groups/OE0436/data/kurilov/nikolas/de_novo_colo800/gffcompare_output/gffCompare.annotated.sorted_ver2.gtf --ref_genome /omics/groups/OE0436/data/kurilov/nikolas/de_novo_rna2/reference_data/hg19.fa --sample A375