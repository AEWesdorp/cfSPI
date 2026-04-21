# Oxford NanoPore pipeline
# Specify snakemake version
from snakemake.utils import min_version
from os.path import join as opj
min_version('9.0.0')

#configfile: 'config/config.yaml'

# This defines the sample name should never contain underscore to avoid confusing rule determination.
wildcard_constraints:
    sample_name='[^_\\W]+'

report: 'report/workflow.rst'

# import functions
include: 'rules/common.smk'

print('Unit file path:', config['units'])

# output directory used to store all results
OUTDIR = config['outdir'] + '/' + config['output_folder']

# Setting `workdir` to the output directory makes sure that the .snakemake caching
# directory is placed here instead of the directory from which the pipeline is 
# launched (usually the git checkout dir), which could fly under the radar and use
# up loads of storage
workdir: config['outdir'] + '/' + config['output_folder']

## sort memory and disk  requirement 
def get_mem_mb(wildcards, attempt):
    return attempt * 16000

def get_disk_mb(wildcards, attempt):
    return attempt * 150000

def get_time(wildcards, attempt):
    return attempt * 60 + 60

def get_time_180_60(wildcards, attempt):
    return attempt * 60 + 180

def get_time_30_120(wildcards, attempt):
    return attempt * 120 + 30


# rules that doesn't require much computational time and power are defied as local rules
localrules: all, get_version_control

# all the output files should be defined at rule all
rule all:
    input:
        # git version
        expand('{OUTDIR}/git-version.log', OUTDIR=OUTDIR),
        # concat_fastq
        expand('{OUTDIR}/results/tmp/{sample_name}.fastq.gz', sample_name=units['sample_name'], OUTDIR=OUTDIR),
        # fastqc
        expand('{OUTDIR}/done/e_fastqc/{sample_name}.done', sample_name=units['sample_name'], OUTDIR=OUTDIR),
        # mapping to host 
        expand('{OUTDIR}/results/host_mapping/{sample_name}_unmapped_host.fastq.gz', sample_name=units['sample_name'],
               OUTDIR=OUTDIR),
        # kraken output
        expand('{OUTDIR}/results/kraken2_report/after_host_mapping/{sample_name}_{database}_conf{k2_threshold}.report',
               database=config['database'], k2_threshold=config['k2_threshold'], sample_name=units['sample_name'],
               OUTDIR=OUTDIR),


rule get_version_control:
    output:
        opj(OUTDIR, 'git-version.log')
    shell:
        """
        echo git branch: > {output};
        git branch >> {output};
        echo ================================ >> {output};
        echo git log: >> {output};
        git log -1  >> {output};
        echo ================================ >> {output};
        echo git status: >> {output};
        git status >> {output};
        echo ================================ >> {output};
        echo git diff: >> {output};
        git diff  >> {output};
        """

rule a_concat_fastq:
    input:
        get_fq
    output:
        fq_raw = temp(opj(OUTDIR, 'results/tmp/{sample_name}.fastq.gz')),
        fq_raw_stats = opj(OUTDIR, 'results/stats/{sample_name}_01_raw_fastq.txt'),
        done = touch(opj(OUTDIR, 'done/a_concat_fastq/{sample_name}.done')),
    benchmark:
        opj(OUTDIR, 'benchmark/a_concat_fastq/{sample_name}.tsv'),
    resources:
        mem_mb  = 8000,
        runtime = get_time,
        disk_mb = 12000,
    threads: 1
    shell:
        """
        mkdir -p {OUTDIR}/results/tmp;
        cat '{input}'/*.fastq.gz > '{output.fq_raw}'     
        echo $(zcat '{output.fq_raw}' | wc -l) / 4 | bc > '{output.fq_raw_stats}';
        """

# rule b_nanoplot:
#     input:
#         fq_raw = rules.a_concat_fastq.output.fq_raw,
#     output:
#         plot_dir = directory(opj(OUTDIR, 'results/stats/nanoplot')),
#         done = touch(opj(OUTDIR, 'done/b_nanoplot/{sample_name}.done')),
#     priority: 47
#     resources:
#         mem_per_cpu = 4000,
#         runtime     = get_time,
#     log:
#         opj(OUTDIR, 'log/b_nanoplot/{sample_name}.log'),
#     threads: 12
#     benchmark:
#         opj(OUTDIR, 'benchmark/b_nanoplot/{sample_name}.tsv'),
#     conda:
#         'envs/nanoplot.yaml'
#     shell:
#         """
#         NanoPlot \
#             --threads {threads} \
#             --fastq {input.fq_raw} \
#             --N50 \
#             --outdir {output.plot_dir}
#         """

rule c_fastplong:
    input:
        fq_raw_c = rules.a_concat_fastq.output.fq_raw,
    output:
        fq_trimmed = temp(opj(OUTDIR, 'results/tmp/{sample_name}.trimmed.fastq.gz')),
        fq_stats = opj(OUTDIR, 'results/stats/{sample_name}_03_fastplong_fastq.txt'),
        html = opj(OUTDIR, 'reports/fastplong/{sample_name}_fastplong.html'),
        json = opj(OUTDIR, 'reports/fastplong/{sample_name}_fastplong.json'),
        failed = temp(opj(OUTDIR, 'results/tmp/fastplong_{sample_name}.fastq.gz')),
        done = touch(opj(OUTDIR, 'done/c_fastplong/{sample_name}.done')),
    priority: 47
    resources:
        mem_mb  = 10000,
        runtime = get_time_30_120,
        disk_mb = 12000,
    log:
        log = opj(OUTDIR, 'log/c_fastplong/{sample_name}.log'),
    threads: 16
    benchmark:
        opj(OUTDIR, 'benchmark/c_fastplong/{sample_name}.tsv'),
    conda:
        'envs/fastplong.yaml'
    shell:
        """
        fastplong \
            --thread {threads} \
            --in {input.fq_raw_c} \
            --out {output.fq_trimmed} \
            --failed_out {output.failed} \
            --html {output.html} \
            --json {output.json} 2>&1 > {log.log};

        echo $(zcat {output.fq_trimmed} | wc -l ) / 4 | bc  > '{output.fq_stats}';
        """

rule e_fastqc:
    input:
        fq_raw_e = rules.a_concat_fastq.output.fq_raw,
        fq_trimmed_e = rules.c_fastplong.output.fq_trimmed,
    output:
        done = touch(opj(OUTDIR, 'done/e_fastqc/{sample_name}.done')),
    priority: 48
    log:
        log = opj(OUTDIR, 'log/e_fastqc/{sample_name}.log'),
    params:
        out_dir = opj(OUTDIR, 'reports/fastqc'),
    resources:
        mem_mb  = 30000,
        runtime = get_time,
        disk_mb = 30000,
    threads: 16
    benchmark:
        opj(OUTDIR, 'benchmark/e_fastqc/{sample_name}.tsv'),
    conda:
        'envs/fastqc.yaml'
    shell:
        """
        mkdir -p {params.out_dir}
        fastqc \
            -t {threads} \
            -o {params.out_dir} \
            {input.fq_raw_e} 2>&1 > {log.log};
        fastqc \
            -t {threads} \
            -o {params.out_dir} \
            {input.fq_trimmed_e} 2>&1 > {log.log};
        """

rule f_host_mapping:
    input:
        fq_trimmed_f = rules.c_fastplong.output.fq_trimmed, 
    output:
        host_bam = temp(opj(OUTDIR, 'results/tmp/{sample_name}_aligned_host.bam')),
        host_unmapped = opj(OUTDIR,'results/host_mapping/{sample_name}_unmapped_host.fastq.gz'),
        fq_stats = opj(OUTDIR,'results/stats/{sample_name}_05_GRCh38_host_mapp_fastq.txt'),
        done = touch(opj(OUTDIR,'done/f_host_mapping/{sample_name}.done')),
    log:
        log = opj(OUTDIR, 'log/f_host_mapping/{sample_name}.log'),
    priority: 49
    params:
        reference = config['reference_genome_dir'] + '/' + config['reference_genome'],
    benchmark:
        opj(OUTDIR, 'benchmark/f_host_mapping/{sample_name}.tsv'),
    conda:
        'envs/minimap2.yaml'
    resources:
        mem_mb_per_cpu  = 2000,
        runtime         = get_time,
        disk_mb         = 50000,
    threads: 32
    shell:
        """
        # Note: minimap2 pipes into samtools view
        minimap2 \
            -x map-ont \
            -t {threads} \
            -a \
            -Y \
            --secondary=yes \
            {params.reference} \
            {input.fq_trimmed_f} | \
        samtools view \
            -b \
            -o {output.host_bam} - ;

        samtools bam2fq \
            --require-flags 4 \
            --threads {threads} \
            -0 {output.host_unmapped} \
            {output.host_bam};
        
        zgrep -c "^@" {output.host_unmapped} > {output.fq_stats};
        """

rule g_kraken2:
    input:
        host_unmapped_g = rules.f_host_mapping.output.host_unmapped,
    output:
        k2_report_hm = opj(OUTDIR, 'results/kraken2_report/after_host_mapping/{sample_name}_{database}_conf{k2_threshold}.report'),
        k2_output_hm = temp(opj(OUTDIR, 'results/kraken2_output/after_host_mapping/{sample_name}_{database}_conf{k2_threshold}.output')),
        k2_output_class_hm = opj(OUTDIR, 'results/kraken2_output/after_host_mapping/{sample_name}_{database}_conf{k2_threshold}.output_classified'),
        done = touch(opj(OUTDIR, 'done/g_kraken2/{sample_name}_{database}_conf{k2_threshold}.done')),
    log:
        log = opj(OUTDIR, 'log/g_kraken2/{sample_name}_{database}_conf{k2_threshold}.log'),
    priority: 50
    params:
        db = config['database_dir'] + '/' + '{database}',
    benchmark:
        opj(OUTDIR, 'benchmark/g_kraken2/{sample_name}_{database}_conf{k2_threshold}.tsv'),
    resources:
        mem_mb  = 128000,
        runtime = 60,
        disk_mb = 128000,
    threads: 16
    conda:
        'envs/kraken.yaml'
    shell:
        """
        kraken2 \
            --confidence {wildcards.k2_threshold} \
            --db {params.db} \
            --threads {threads} \
            --report-zero-counts \
            --report-minimizer-data \
            --output {output.k2_output_hm} \
            --report {output.k2_report_hm} \
            {input.host_unmapped_g};

        grep -v "^U" {output.k2_output_hm} > {output.k2_output_class_hm};
        """