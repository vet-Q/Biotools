import subprocess


def bcftools_view(input_vcf, output_vcf, dry_run=False, **kwargs):
    """
    Run `bcftools view` on an `input_vcf`

    params:
        input_vcf: str
            Path to input `.vcf`.
        output_vcf: str
            Path to output `.vcf`.
        dry_run: bool
            Print command instead of running it.
        kwargs: key=value
            Additional arguments will be passed to
            bcftools view as flags; e.g. 
            `-<key> <value>`. 
        
    returns
        None

    """

    cmd = f"bcftools view -I {input_vcf} "
    cmd += " ".join([f"-{k} {v}" for k, v in kwargs.items()])
    cmd += f" -o {output_vcf}"

    if dry_run:
        print(cmd)
        return

    subprocess.run(cmd, shell=True, check=True)


def bcftools_index(input_vcf):
    """
    Run `bcftools index`
    
    params
        input_vcf : str
            VCF file to be indexed.
    
    returns
        None

    """
    
    cmd = f"bcftools index {input_vcf}" 
    subprocess.run(cmd, shell=True, check=True)


def bcftools_merge(input_vcfs, output_vcf, dry_run=False):
    """
    Run `bcftools merge`
    
    params
        input_vcfs: list[str]
            List of paths to `.vcf` files that will be merged.
            Note that these need to be compressed for merge
            to work.
        output_vcf: str
            Path to output vcf.
        dry_run: bool
            Print command instead of running.
    
    """
    cmd = "bcftools merge"
    cmd += f" {' '.join(input_vcfs)}"
    cmd += f" -o {output_vcf}"

    if dry_run:
        print(cmd)
        return

    subprocess.run(cmd, shell=True, check=True)
