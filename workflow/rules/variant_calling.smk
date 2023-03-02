rule pull_docker_VC:
    output:
        touch(resolve_results_filepath(
            config.get("paths").get("results_dir"),"ctat_container.pull.done"))
    params:
        container=config.get("docker").get("ctat")
    message: "Downloading the following Docker container: {params.container}."
    shell:
        "docker pull {params.container} "


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
        container=config.get("docker").get("ctat"),
        results_dir=config.get("paths").get("results_dir"),
        ctat_path=config.get("resources").get("ctat_path"),
        fq_path=resolve_results_filepath(
            config.get("paths").get("results_dir"),"reads/trimmed/"),
        r1_name="{sample}-R1-trimmed.fq.gz",
        r2_name="{sample}-R2-trimmed.fq.gz",
        outprefix="variant_calling/{sample}",
        sample_id="{sample}"
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),"logs/docker_vc/{sample}.log")
    threads: conservative_cpu_count(reserve_cores=1, max_cores=10)
    resources:
        tmpdir=config.get("paths").get("tmp_dir")
    message: "Running Docker container for Variant Calling with {threads} threads for the following files {input.read1}{input.read2}."
    shell:
        "docker run -v {params.results_dir}:/data "
        "-v {params.ctat_path}:/ctat_genome_lib_build_dir "
        "-v {params.fq_path}:/reads --rm "
        "{params.container} "
        "/usr/local/src/ctat-mutations/ctat_mutations "
        "--left /reads/{params.r1_name} "
        "--right /reads/{params.r2_name} "
        "--genome_lib_dir /ctat_genome_lib_build_dir "
        "--outputdir /data/{params.outprefix} "
        "--sample_id {params.sample_id} "
        "--CPU {threads} "
        ">& {log} "