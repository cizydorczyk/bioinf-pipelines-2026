import os
import sys

##############################################
# Configuration validation and setup
##############################################

# Required config parameters
REQUIRED_CONFIG = ['fastq_dir', 'project_dir', 'ref', 'st', 'isolate_list', 
                   'fastq_endings', 'remove_ref', 'outgroups', 'snippy_minfrac',
                   'snippy_cleanup', 'snippy_unmapped', 'model', 'bb', 'mask_file']

# Validate required config parameters
missing_params = [param for param in REQUIRED_CONFIG if param not in config]
if missing_params:
    raise ValueError(f"Missing required config parameters: {', '.join(missing_params)}")

# Extract config with defaults
fastq_dir = config['fastq_dir']
project_dir = config['project_dir']
ref = config['ref']
st = config['st']
isolate_list_handle = config['isolate_list']
fastq_endings = config['fastq_endings'].split(",")
remove_ref = config['remove_ref']
outgroups = config['outgroups']

snippy_minfrac = config['snippy_minfrac']
snippy_cleanup = config['snippy_cleanup']
snippy_unmapped = config['snippy_unmapped']

iqtree_model = config['model']
iqtree_bb = config['bb']

mask_file = config['mask_file']

# Threading configuration
# total_threads: total threads available to pipeline (e.g., 24)
# threads_per_snippy: threads per snippy job (e.g., 4)
# This allows: 24/4 = 6 snippy jobs to run simultaneously
total_threads = config.get('threads', 24)
threads_per_snippy = config.get('threads_per_snippy', 4)

# Validate threading configuration
# if total_threads % threads_per_snippy != 0:
#     raise ValueError(
#         f"threads_per_snippy ({threads_per_snippy}) must evenly divide "
#         f"total threads ({total_threads}). "
#         f"Valid values: {[i for i in range(1, total_threads+1) if total_threads % i == 0]}"
#     )
if total_threads % threads_per_snippy != 0:
    valid_values = [i for i in range(1, total_threads+1) if total_threads % i == 0]
    raise ValueError(
        f"threads_per_snippy ({threads_per_snippy}) must evenly divide "
        f"total threads ({total_threads}). "
        f"Valid values: {valid_values}"
    )

max_parallel_snippy = total_threads // threads_per_snippy

# Validate inputs exist
if not os.path.exists(fastq_dir):
    raise FileNotFoundError(f"FASTQ directory not found: {fastq_dir}")
if not os.path.exists(ref):
    raise FileNotFoundError(f"Reference file not found: {ref}")
if not os.path.exists(isolate_list_handle):
    raise FileNotFoundError(f"Isolate list not found: {isolate_list_handle}")

##############################################
# Global wildcard constraints
##############################################

wildcard_constraints:
    sample = '[a-zA-Z0-9_:;-]+',
    st = '[a-zA-Z0-9_:;]+',
    remove_seqs_flag = '[a-zA-Z0-9_.]+',
    remove_seqs_flag_no_dot = '[a-zA-Z0-9_]+'

##############################################
# Helper functions
##############################################

def determine_removal_flags():
    """
    Determine which sequence removal flags to use based on config.
    Returns: tuple of (flags_with_dot, flags_without_dot)
    
    CRITICAL: This must produce IDENTICAL output to original logic.
    """
    has_remove_ref = remove_ref.split(",")[0] != ""
    has_outgroups = outgroups.split(",")[0] != ""
    
    if not has_remove_ref and not has_outgroups:
        # Only run with reference
        return ([""], [""])
    elif has_remove_ref and not has_outgroups:
        # Remove reference only
        return (["no_ref."], ["no_ref"])
    elif not has_remove_ref and has_outgroups:
        # Run with and without outgroups
        return (["with_outgroups.", "without_outgroups."], 
                ["with_outgroups", "without_outgroups"])
    else:
        # Remove reference and run with/without outgroups
        return (["no_ref_with_outgroups.", "no_ref_without_outgroups."],
                ["no_ref_with_outgroups", "no_ref_without_outgroups"])

def load_isolates():
    """Load isolate list from file."""
    with open(isolate_list_handle, 'r') as f:
        isolates = [line.strip() for line in f if line.strip()]
    
    if not isolates:
        raise ValueError(f"No isolates found in {isolate_list_handle}")
    
    return isolates

def get_isolates_with_outgroups(isolates):
    """
    Combine isolates with outgroups if specified.
    CRITICAL: Must match original logic exactly.
    """
    if outgroups != "":
        return isolates + outgroups.split(",")
    return isolates

def build_fastq_dict(isolates_all):
    """
    Build dictionary mapping isolates to their FASTQ files.
    CRITICAL: Must produce same dict structure as original.
    Validates that all required files exist.
    """
    fastq_files = [os.path.abspath(os.path.join(fastq_dir, fq)) for fq in os.listdir(fastq_dir)]
    fastq_dict = {}
    
    for isolate in isolates_all:
        # Find matching forward and reverse reads
        forward = [fastq for fastq in fastq_files if isolate in fastq and fastq_endings[0] in fastq]
        reverse = [fastq for fastq in fastq_files if isolate in fastq and fastq_endings[1] in fastq]
        
        # Validate files were found
        if not forward:
            raise FileNotFoundError(
                f"No forward reads found for isolate '{isolate}' "
                f"with ending '{fastq_endings[0]}' in {fastq_dir}"
            )
        if not reverse:
            raise FileNotFoundError(
                f"No reverse reads found for isolate '{isolate}' "
                f"with ending '{fastq_endings[1]}' in {fastq_dir}"
            )
        if len(forward) > 1:
            raise ValueError(
                f"Multiple forward read files found for isolate '{isolate}': {forward}"
            )
        if len(reverse) > 1:
            raise ValueError(
                f"Multiple reverse read files found for isolate '{isolate}': {reverse}"
            )
        
        # CRITICAL: Use same dict structure as original
        fastq_dict[isolate] = {'R1': forward[0], 'R2': reverse[0]}
    
    return fastq_dict

##############################################
# Initialize variables - MUST match original
##############################################

remove_seqs_flag, remove_seqs_flag_no_dot = determine_removal_flags()
isolate_list = load_isolates()
isolate_list_with_outgroups = get_isolates_with_outgroups(isolate_list)
fastq_files_dict = build_fastq_dict(isolate_list_with_outgroups)

# Mask file parameter - MUST match original
if len(mask_file) > 1:
    maskfile = f" --mask {mask_file}"
else:
    maskfile = ""

# Report configuration
print(f"Pipeline configuration:", file=sys.stderr)
print(f"  - Isolates: {len(isolate_list)}", file=sys.stderr)
print(f"  - With outgroups: {len(isolate_list_with_outgroups)}", file=sys.stderr)
print(f"  - Removal flags: {remove_seqs_flag_no_dot}", file=sys.stderr)
print(f"  - Total threads: {total_threads}", file=sys.stderr)
print(f"  - Threads per snippy job: {threads_per_snippy}", file=sys.stderr)
print(f"  - Max parallel snippy jobs: {max_parallel_snippy}", file=sys.stderr)

# CRITICAL: Create directories exactly as original does
# (Snakemake doesn't auto-create all intermediate dirs reliably)
if not os.path.isdir(project_dir):
    os.mkdir(project_dir)

if not os.path.isdir(os.path.join(project_dir, "snippy-core")):
    os.mkdir(os.path.join(project_dir, "snippy-core"))

if not os.path.isdir(os.path.join(project_dir, "iqtree")):
    os.mkdir(os.path.join(project_dir, "iqtree"))

if not os.path.isdir(os.path.join(project_dir, "cfml")):
    os.mkdir(os.path.join(project_dir, "cfml"))

if not os.path.isdir(os.path.join(project_dir, "logs")):
    os.mkdir(os.path.join(project_dir, "logs"))

##############################################
# Input functions
##############################################

def get_reads(wildcards):
    """Get paired-end reads for a sample."""
    return {
        'R1': fastq_files_dict[wildcards.sample]["R1"],
        'R2': fastq_files_dict[wildcards.sample]["R2"]
    }

##############################################
# Target rules - MUST produce identical outputs
##############################################

rule all:
    input:
        #### snippy: ####
        expand("{output}raw-snippy-output/{sample}", output=project_dir, sample=isolate_list_with_outgroups),
        #### snippy-core: ####
        expand("{output}snippy-core/{st}.full.aln", output=project_dir, st=st),
        #### snippy-clean_full_aln: ####
        expand("{output}snippy-core/{st}.clean.full.aln", output=project_dir, st=st),
        #### remove_seq: ####
        expand("{output}snippy-core/{st}.{remove_seqs_flag}clean.full.aln", output=project_dir, st=st, remove_seqs_flag=remove_seqs_flag),
        #### iqtree - only for >=3 isolates ####
        expand("{output}iqtree/{st}.{remove_seqs_flag_no_dot}.contree", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot) if len(isolate_list) >= 4 else (expand("{output}iqtree/{st}.{remove_seqs_flag_no_dot}.treefile", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot) if len(isolate_list) == 3 else []),
        #### cfml - only for >=3 isolates ####
        expand("{output}cfml/{st}.{remove_seqs_flag_no_dot}.labelled_tree.newick", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot) if len(isolate_list) >= 3 else [],
        #### maskrcsvg - only for >=3 isolates ####
        expand("{output}cfml/{st}.rc_masked.{remove_seqs_flag_no_dot}.clean.full.aln", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot) if len(isolate_list) >= 3 else [],
        #### snpdists - conditional based on isolate count ####
        expand("{output}cfml/{st}.direct.{remove_seqs_flag_no_dot}.clean.full.aln.snpdists_matrix.txt", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot) if len(isolate_list) == 2 else expand("{output}cfml/{st}.rc_masked.{remove_seqs_flag_no_dot}.clean.full.aln.snpdists_matrix.txt", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot)

# Individual target rules for partial pipeline execution
rule run_snippy:
    input:
        expand("{output}raw-snippy-output/{sample}",output=project_dir,sample=isolate_list_with_outgroups)

rule run_snippy_core:
    input:
        expand("{output}snippy-core/{st}.full.aln",output=project_dir,st=st)

rule run_snippy_clean:
    input:
        expand("{output}snippy-core/{st}.clean.full.aln",output=project_dir,st=st)

rule run_remove_seqs:
    input:
        expand("{output}snippy-core/{st}.{remove_seqs_flag}clean.full.aln", output=project_dir, st=st, remove_seqs_flag=remove_seqs_flag)

rule run_iqtree:
    input:
        expand("{output}iqtree/{st}.{remove_seqs_flag_no_dot}.contree",output=project_dir,st=st,remove_seqs_flag_no_dot=remove_seqs_flag_no_dot) if len(isolate_list) >= 4 else expand("{output}iqtree/{st}.{remove_seqs_flag_no_dot}.treefile",output=project_dir,st=st,remove_seqs_flag_no_dot=remove_seqs_flag_no_dot)

rule run_cfml:
    input:
        expand("{output}cfml/{st}.{remove_seqs_flag_no_dot}.labelled_tree.newick",output=project_dir,st=st,remove_seqs_flag_no_dot=remove_seqs_flag_no_dot)

rule run_maskrcsvg:
    input:
        expand("{output}cfml/{st}.rc_masked.{remove_seqs_flag_no_dot}.clean.full.aln",output=project_dir,st=st,remove_seqs_flag_no_dot=remove_seqs_flag_no_dot)

rule run_snpdists:
    input:
        expand("{output}cfml/{st}.direct.{remove_seqs_flag_no_dot}.clean.full.aln.snpdists_matrix.txt", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot) if len(isolate_list) == 2 else expand("{output}cfml/{st}.rc_masked.{remove_seqs_flag_no_dot}.clean.full.aln.snpdists_matrix.txt", output=project_dir, st=st, remove_seqs_flag_no_dot=remove_seqs_flag_no_dot)

##############################################
# Pipeline rules
##############################################

rule snippy:
    """
    Run Snippy variant calling on paired-end reads.
    
    Threading: Uses threads_per_snippy threads per job.
    Snakemake will run max_parallel_snippy jobs simultaneously.
    Example: 24 total threads, 4 per job = 6 parallel jobs
    """
    input:
        unpack(get_reads)
    output:
        directory("{output}raw-snippy-output/{sample}")
    params:
        reference=ref,
        minfrac=snippy_minfrac,
        cleanup = " --cleanup" if snippy_cleanup == "yes" else "",
        unmapped = " --unmapped" if snippy_unmapped == "yes" else ""
    threads: threads_per_snippy
    log:
        "{output}logs/snippy/{sample}.log"
    shell:
        """
        snippy \
            --R1 {input.R1} \
            --R2 {input.R2} \
            --ref {params.reference} \
            --outdir {output} \
            --cpus {threads} \
            --minfrac {params.minfrac} \
            {params.cleanup} \
            {params.unmapped} \
            2>&1 | tee {log}
        """

rule snippy_core:
    """
    Generate core SNP alignment from Snippy outputs.
    CRITICAL: Output structure must match original exactly.
    """
    input:
        expand("{output}raw-snippy-output/{sample}", output=project_dir, sample=isolate_list_with_outgroups)
    output:
        "{output}snippy-core/{st}.full.aln"
    params:
        reference=ref,
        output_prefix="{output}snippy-core/{st}",
        mf=maskfile
    log:
        "{output}logs/snippy-core/{st}.log"
    shell:
        """
        snippy-core \
            --prefix {params.output_prefix} \
            {params.mf} \
            --ref {params.reference} \
            {input} \
            2>&1 | tee {log}
        """

rule snippy_clean:
    """
    Clean the full alignment to remove non-variant sites.
    """
    input:
        "{output}snippy-core/{st}.full.aln"
    output:
        "{output}snippy-core/{st}.clean.full.aln"
    log:
        "{output}logs/snippy-clean/{st}.log"
    shell:
        """
        snippy-clean_full_aln {input} > {output} 2> {log}
        """

rule remove_seqs:
    """
    Remove reference and/or outgroup sequences from alignment.
    CRITICAL: Must produce EXACT same outputs as original.
    
    Logic matches original remove_seqs_function.py exactly:
    - No remove_ref, no outgroups: Keep everything (empty flag)
    - Remove_ref only: Remove ref (no_ref flag)
    - Outgroups only: Create with/without outgroups variants
    - Both: Remove ref, create with/without outgroups variants
    """
    input:
        expand("{output}snippy-core/{st}.clean.full.aln", output=project_dir, st=st)
    output:
        expand("{output}snippy-core/{st}.{remove_seqs_flag}clean.full.aln", output=project_dir, st=st, remove_seqs_flag=remove_seqs_flag)
    params:
        remove_ref=remove_ref,
        remove_outgroups=outgroups
    log:
        expand("{output}logs/remove-seqs/{st}.log", output=project_dir, st=st)
    run:
        from Bio import SeqIO
        import sys
        
        def remove_seqs(seq_to_remove, fasta, output_fasta):
            """Remove specified sequences from FASTA file."""
            records = list(SeqIO.parse(fasta, "fasta"))
            records2 = [i for i in records if i.id not in seq_to_remove.split(",")]
            
            with open(output_fasta, "w") as outfile:
                SeqIO.write(records2, outfile, "fasta")
        
        # Control flow for removing ref/outgroups (EXACT COPY of original logic):
        if params.remove_ref.split(",")[0] == "" and params.remove_outgroups.split(",")[0] == "":
            # Only run with ref if neither remove ref nor outgroups set
            print("Not removing anything, running with clean.{st}.full.aln...")
            remove_seqs("", input[0], output[0])
        
        elif params.remove_ref.split(",")[0] != "" and params.remove_outgroups.split(",")[0] == "":
            # If remove ref set, run without ref (not with)
            print("Removing ref and running without ref only...")
            remove_seqs(params.remove_ref, input[0], output[0])
        
        elif params.remove_ref.split(",")[0] == "" and params.remove_outgroups.split(",")[0] != "":
            # If outgroups present, always run with and without them
            print("Keeping ref and running with/without outgroups...")
            remove_seqs("", input[0], output[0])
            remove_seqs(params.remove_outgroups, input[0], output[1])
        
        elif params.remove_ref.split(",")[0] != "" and params.remove_outgroups.split(",")[0] != "":
            # If outgroups present and remove ref, remove ref & run with and without outgroups
            print("Removing ref and running with/without outgroups...")
            remove_seqs(params.remove_ref, input[0], output[0])
            remove_seqs(params.remove_ref + "," + params.remove_outgroups, input[0], output[1])
        
        else:
            print("Check remove reference/outgroups flags and try again.")
            sys.exit()

rule iqtree:
    """
    Infer maximum likelihood phylogeny with IQ-TREE.
    CRITICAL: Output extension changes based on number of isolates.
    Only runs for >=3 isolates.
    """
    input:
        "{output}snippy-core/{st}.{remove_seqs_flag_no_dot}.clean.full.aln"
    output:
        "{output}iqtree/{st}.{remove_seqs_flag_no_dot}.contree" if len(isolate_list) >= 4 else "{output}iqtree/{st}.{remove_seqs_flag_no_dot}.treefile"
    params:
        iqtree_params = f"-m {iqtree_model} -bb {iqtree_bb}" if len(isolate_list) >= 4 else f"-m {iqtree_model}",
        iqtree_prefix = "{output}iqtree/{st}.{remove_seqs_flag_no_dot}"
    threads: total_threads
    log:
        "{output}logs/iqtree/{st}.{remove_seqs_flag_no_dot}.log"
    shell:
        """
        iqtree \
            -s {input} \
            {params.iqtree_params} \
            -pre {params.iqtree_prefix} \
            -nt {threads} \
            2>&1 | tee {log}
        """

rule clonalframeml:
    """
    Detect recombination with ClonalFrameML.
    Only runs for >=3 isolates.
    """
    input:
        tree = "{output}iqtree/{st}.{remove_seqs_flag_no_dot}.contree" if len(isolate_list) >= 4 else "{output}iqtree/{st}.{remove_seqs_flag_no_dot}.treefile",
        fasta = "{output}snippy-core/{st}.{remove_seqs_flag_no_dot}.clean.full.aln"
    output:
        "{output}cfml/{st}.{remove_seqs_flag_no_dot}.labelled_tree.newick"
    params:
        cfml_prefix = "{output}cfml/{st}.{remove_seqs_flag_no_dot}"
    log:
        "{output}logs/cfml/{st}.{remove_seqs_flag_no_dot}.log"
    shell:
        """
        ClonalFrameML \
            {input.tree} \
            {input.fasta} \
            {params.cfml_prefix} \
            -em true \
            -show_progress true \
            2>&1 | tee {log}
        """

rule maskrcsvg:
    """
    Mask recombinant regions from alignment.
    CRITICAL: cfml_tree input forces execution order (after cfml completes).
    Only runs for >=3 isolates.
    """
    input:
        fasta = "{output}snippy-core/{st}.{remove_seqs_flag_no_dot}.clean.full.aln",
        cfml_tree = "{output}cfml/{st}.{remove_seqs_flag_no_dot}.labelled_tree.newick"
    output:
        "{output}cfml/{st}.rc_masked.{remove_seqs_flag_no_dot}.clean.full.aln"
    params:
        cfml_output_prefix = "{output}cfml/{st}.{remove_seqs_flag_no_dot}",
        maskrc_regions_file = "{output}cfml/{st}.rc_masked.{remove_seqs_flag_no_dot}.regions.txt"
    log:
        "{output}logs/maskrcsvg/{st}.{remove_seqs_flag_no_dot}.log"
    shell:
        """
        maskrc-svg.py \
            --symbol N \
            --aln {input.fasta} \
            --out {output} \
            --regions {params.maskrc_regions_file} \
            {params.cfml_output_prefix} \
            2>&1 | tee {log}
        """

rule snpdists_direct:
    """
    Calculate pairwise SNP distance matrix directly from cleaned alignment.
    Used when isolate count is too low for phylogenetic analysis (2 isolates).
    """
    input:
        "{output}snippy-core/{st}.{remove_seqs_flag_no_dot}.clean.full.aln"
    output:
        "{output}cfml/{st}.direct.{remove_seqs_flag_no_dot}.clean.full.aln.snpdists_matrix.txt"
    log:
        "{output}logs/snpdists-direct/{st}.{remove_seqs_flag_no_dot}.log"
    shell:
        """
        snp-dists {input} > {output} 2> {log}
        """

rule snpdists:
    """
    Calculate pairwise SNP distance matrix from recombination-masked alignment.
    Used when isolate count allows phylogenetic analysis (>=3 isolates).
    """
    input:
        "{output}cfml/{st}.rc_masked.{remove_seqs_flag_no_dot}.clean.full.aln"
    output:
        "{output}cfml/{st}.rc_masked.{remove_seqs_flag_no_dot}.clean.full.aln.snpdists_matrix.txt"
    log:
        "{output}logs/snpdists/{st}.{remove_seqs_flag_no_dot}.log"
    shell:
        """
        snp-dists {input} > {output} 2> {log}
        """