#!/usr/bin/env python
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#Initial-setup" data-toc-modified-id="Initial-setup-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Initial setup</a></span></li><li><span><a href="#dsub-test" data-toc-modified-id="dsub-test-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>dsub test</a></span><ul class="toc-item"><li><span><a href="#Job-script" data-toc-modified-id="Job-script-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Job script</a></span></li><li><span><a href="#Run-dsub" data-toc-modified-id="Run-dsub-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>Run dsub</a></span></li><li><span><a href="#10-sample-batch-runtime" data-toc-modified-id="10-sample-batch-runtime-2.3"><span class="toc-item-num">2.3&nbsp;&nbsp;</span>10-sample batch runtime</a></span></li><li><span><a href="#3-sample-batch-runtime" data-toc-modified-id="3-sample-batch-runtime-2.4"><span class="toc-item-num">2.4&nbsp;&nbsp;</span>3-sample batch runtime</a></span></li><li><span><a href="#1-sample-batch-runtime" data-toc-modified-id="1-sample-batch-runtime-2.5"><span class="toc-item-num">2.5&nbsp;&nbsp;</span>1-sample batch runtime</a></span></li><li><span><a href="#kir-dsub-outputs" data-toc-modified-id="kir-dsub-outputs-2.6"><span class="toc-item-num">2.6&nbsp;&nbsp;</span>kir dsub outputs</a></span></li></ul></li></ul></div>

# ## Initial setup

# In[ ]:


get_ipython().system('pip install tctk --upgrade')


# In[ ]:


get_ipython().system('pip show tctk')


# In[ ]:


from tctk import AoUTools as at
import os


# In[ ]:


import os
import polars as pl


# In[ ]:


my_bucket = os.getenv('WORKSPACE_BUCKET')
ancestry_df = pl.read_csv(f'{my_bucket}/data/ancestry_metadata.tsv', separator='\t')


# In[ ]:


ancestry_df.group_by('ancestry_pred_other').len()


# In[ ]:


bucket = os.getenv("WORKSPACE_BUCKET")
bucket


# In[ ]:


get_ipython().system('gsutil -u $GOOGLE_PROJECT ls gs://fc-aou-datasets-controlled/v8/wgs/cram/')


# In[ ]:


get_ipython().system('gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v8/wgs/cram/manifest.csv .')


# In[ ]:


get_ipython().system('head -n 5 manifest.csv')


# In[ ]:


get_ipython().system('gsutil -u $GOOGLE_PROJECT ls gs://fc-aou-datasets-controlled/pooled/wgs/cram/v8_delta/wgs_1000000*')


# ## dsub test

# ### Job script

# In[ ]:


get_ipython().run_cell_magic('writefile', 'kir_dsub_test.sh', '\n#!/bin/bash\n\n# select a region and convert cram to bam\n# this is without looping, using 3 participants to test scaling up to 10/VM; try looping for large scale\n# 01\ngatk PrintReads \\\n-R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\\n-I $INPUT_CRAM_1 \\\n--read-index $INPUT_CRAI_1 \\\n-L chr19:54000000-55100000 \\\n--gcs-project-for-requester-pays $GOOGLE_PROJECT \\\n-O sample01_chr19_region.bam\n        \n# 02\ngatk PrintReads \\\n-R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\\n-I $INPUT_CRAM_2 \\\n--read-index $INPUT_CRAI_2 \\\n-L chr19:54000000-55100000 \\\n--gcs-project-for-requester-pays $GOOGLE_PROJECT \\\n-O sample02_chr19_region.bam\n        \n# 03\ngatk PrintReads \\\n-R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\\n-I $INPUT_CRAM_3 \\\n--read-index $INPUT_CRAI_3 \\\n-L chr19:54000000-55100000 \\\n--gcs-project-for-requester-pays $GOOGLE_PROJECT \\\n-O sample03_chr19_region.bam\n        \n# 04\ngatk PrintReads \\\n-R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\\n-I $INPUT_CRAM_1 \\\n--read-index $INPUT_CRAI_1 \\\n-L chr19:54000000-55100000 \\\n--gcs-project-for-requester-pays $GOOGLE_PROJECT \\\n-O sample04_chr19_region.bam\n        \n# 05\ngatk PrintReads \\\n-R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\\n-I $INPUT_CRAM_2 \\\n--read-index $INPUT_CRAI_2 \\\n-L chr19:54000000-55100000 \\\n--gcs-project-for-requester-pays $GOOGLE_PROJECT \\\n-O sample05_chr19_region.bam\n        \n# 06\ngatk PrintReads \\\n-R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\\n-I $INPUT_CRAM_3 \\\n--read-index $INPUT_CRAI_3 \\\n-L chr19:54000000-55100000 \\\n--gcs-project-for-requester-pays $GOOGLE_PROJECT \\\n-O sample06_chr19_region.bam\n        \n# 07\ngatk PrintReads \\\n-R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\\n-I $INPUT_CRAM_1 \\\n--read-index $INPUT_CRAI_1 \\\n-L chr19:54000000-55100000 \\\n--gcs-project-for-requester-pays $GOOGLE_PROJECT \\\n-O sample07_chr19_region.bam\n        \n# 08\ngatk PrintReads \\\n-R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\\n-I $INPUT_CRAM_2 \\\n--read-index $INPUT_CRAI_2 \\\n-L chr19:54000000-55100000 \\\n--gcs-project-for-requester-pays $GOOGLE_PROJECT \\\n-O sample08_chr19_region.bam\n        \n# 09\ngatk PrintReads \\\n-R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\\n-I $INPUT_CRAM_3 \\\n--read-index $INPUT_CRAI_3 \\\n-L chr19:54000000-55100000 \\\n--gcs-project-for-requester-pays $GOOGLE_PROJECT \\\n-O sample09_chr19_region.bam\n        \n# 10\ngatk PrintReads \\\n-R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\\n-I $INPUT_CRAM_1 \\\n--read-index $INPUT_CRAI_1 \\\n-L chr19:54000000-55100000 \\\n--gcs-project-for-requester-pays $GOOGLE_PROJECT \\\n-O sample10_chr19_region.bam\n        \n            \n# run kir-mapper\n# Step 1: Map (using your BAM file)\nkir-mapper map -bam sample01_chr19_region.bam -sample sample01 -output ./kir_output -threads 3\nkir-mapper map -bam sample02_chr19_region.bam -sample sample02 -output ./kir_output -threads 3\nkir-mapper map -bam sample03_chr19_region.bam -sample sample03 -output ./kir_output -threads 3\nkir-mapper map -bam sample04_chr19_region.bam -sample sample04 -output ./kir_output -threads 3\nkir-mapper map -bam sample05_chr19_region.bam -sample sample05 -output ./kir_output -threads 3\nkir-mapper map -bam sample06_chr19_region.bam -sample sample06 -output ./kir_output -threads 3\nkir-mapper map -bam sample07_chr19_region.bam -sample sample07 -output ./kir_output -threads 3\nkir-mapper map -bam sample08_chr19_region.bam -sample sample08 -output ./kir_output -threads 3\nkir-mapper map -bam sample09_chr19_region.bam -sample sample09 -output ./kir_output -threads 3\nkir-mapper map -bam sample10_chr19_region.bam -sample sample10 -output ./kir_output -threads 3\n\n# subsequent steps can be done with 1 single command for all samples\n# Step 2: Copy number estimation\nkir-mapper ncopy -output ./kir_output -threads 3\n\n# Step 3: Genotype (SNP and allele calling)\nkir-mapper genotype -output ./kir_output -threads 3\n\n# Step 4: Haplotype (requires ≥50 samples for accuracy)\nkir-mapper haplotype -output ./kir_output -threads 3\n\n# Copy kir output folder to dsub output folder\ncp -r kir_output/* $OUTPUT_DIR\n')


# ### Run dsub

# In[ ]:


kir_dsub = at.Dsub(
    provider="google-batch",
    machine_type="c2d-highcpu-4",
    docker_image="phetk/gatk-kirmapper:0.1",
    job_script_name="kir_dsub_test.sh",
    job_name="kir_dsub",
    env_dict={
        "GOOGLE_PROJECT": os.getenv("GOOGLE_PROJECT"),
        "INPUT_CRAM_1": "gs://fc-aou-datasets-controlled/pooled/wgs/cram/v8_delta/wgs_1000000.cram",
        "INPUT_CRAM_2": "gs://fc-aou-datasets-controlled/pooled/wgs/cram/v6_base/wgs_1000004.cram",
        "INPUT_CRAM_3": "gs://fc-aou-datasets-controlled/pooled/wgs/cram/v6_base/wgs_1000033.cram"
    },
    input_dict={
        "INPUT_CRAI_1": "gs://fc-aou-datasets-controlled/pooled/wgs/cram/v8_delta/wgs_1000000.cram.crai",
        "INPUT_CRAI_2": "gs://fc-aou-datasets-controlled/pooled/wgs/cram/v6_base/wgs_1000004.cram.crai",
        "INPUT_CRAI_3": "gs://fc-aou-datasets-controlled/pooled/wgs/cram/v6_base/wgs_1000033.cram.crai"
    },
    output_dict={},
    custom_args=f"--output-recursive OUTPUT_DIR={bucket}/dsub/results/kir/"
)
kir_dsub.run(show_command=True)


# ### 10-sample batch runtime

# In[ ]:


# 10 samples for the record - dont run again
kir_dsub.check_status(streaming=True)


# ### 3-sample batch runtime

# In[ ]:


# 3 samples for the record - dont run again
kir_dsub.check_status(streaming=True)


# ### 1-sample batch runtime

# In[ ]:


# 1 sample for the record - dont run again
kir_dsub.check_status(streaming=True)


# ### kir dsub outputs

# In[ ]:


get_ipython().system('gsutil ls {bucket}/dsub/results/kir/')


# In[ ]:


get_ipython().system('gsutil ls {bucket}/dsub/results/kir/map/')


# In[ ]:


get_ipython().system('gsutil ls {bucket}/dsub/results/kir/ncopy/')


# In[ ]:


get_ipython().system('gsutil cat {bucket}/dsub/results/kir/ncopy/copy_numbers.table.txt | head -n 5')


# In[ ]:


get_ipython().system('gsutil cat {bucket}/dsub/results/kir/ncopy/presence.table.txt | head -n 5')


# In[ ]:


get_ipython().system('gsutil cp {bucket}/dsub/results/kir/ncopy/plots/KIR2DL3.png .')


# In[ ]:


kir_dsub.kill_all()


# In[ ]:




