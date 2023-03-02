import os
rule pull_docker_VC:
    output:
        touch(resolve_results_filepath(
            config.get("paths").get("results_dir"),"ctat_container.pull.done"))
    params:
        container=config.get("docker").get("ctat").get("container")
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/docker/pull_docker.log")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/bash.yaml")
    message: "Downloading the following Docker container: {params.container}."
    shell:
        "docker pull {params.container} "
        ">& {log} "


rule docker_VC:
    input:
        read1=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/{sample}-R1-trimmed.fq.gz"),
        read2=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/{sample}-R2-trimmed.fq.gz"),
        docker=resolve_results_filepath(
            config.get("paths").get("results_dir"),"ctat_container.pull.done")
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"variant_calling/{sample}/{sample}.vcf.gz")
    params:
        container=config.get("docker").get("ctat").get("container"),
        path=config.get("docker").get("ctat").get("path"),
        results_dir=lambda wildcards, output: return_res_dir(output, elements=3),
        ctat_path=config.get("resources").get("ctat_path"),
        fq_path=lambda wildcards, input: os.path.dirname(input.read1),
        r1_name=lambda wildcards, input: os.path.basename(input.read2),
        r2_name=lambda wildcards, input: os.path.basename(input.read1),
        outprefix="variant_calling/{sample}",
        sample_id="{sample}"
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/docker_vc/{sample}.log")
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"),"workflow/envs/bash.yaml")
    threads: conservative_cpu_count(reserve_cores=1, max_cores=10)
    resources:
        tmpdir=config.get("paths").get("tmp_dir")
    message: "Running Docker container for Variant Calling with {threads} threads for the following files {input.read1}{input.read2}."
    shell:
        "docker run -v {params.results_dir}:/data "
        "-v {params.ctat_path}:/ctat_genome_lib_build_dir "
        "-v {params.fq_path}:/reads --rm "
        "{params.container} "
        "{params.path} "
        "--left /reads/{params.r1_name} "
        "--right /reads/{params.r2_name} "
        "--genome_lib_dir /ctat_genome_lib_build_dir "
        "--outputdir /data/{params.outprefix} "
        "--sample_id {params.sample_id} "
        "--CPU {threads} "
        ">& {log} "