# RNA-seq preprocessing

This pipeline processes six public RNA-seq libraries with FastQC/MultiQC, STAR
two-pass alignment, reverse-stranded featureCounts, RegTools junction extraction,
and RSeQC strand assessment. Fresh PE and SE repeats and a featureCounts strand
benchmark reproduce the validation experiments in the accompanying article.

## Run

Requires Windows x86-64, PowerShell 7, CPython 3.13, curl, and Docker Desktop
running Linux containers. The default configuration uses eight threads, a
26-GiB STAR limit, and approximately 32 GiB physical RAM. Allow 40 GiB free for
sample processing, plus space for the shared reference/index and retained results.
Internet access is required to acquire reads, references, containers, and Python
dependencies.

Run from the repository root:

```powershell
python scripts/verify_checksums.py
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --require-hashes -r requirements.lock
pwsh -NoProfile -File scripts/run_pipeline.ps1 -HostData D:/RNAseq-data
```

The complete run prepares the reference, processes six libraries serially,
runs both fresh repeats and the strand benchmark, and generates tables, figures,
MultiQC, resource measurements, and execution manifests in `outputs/`.
Reads and intermediates use `-HostData`; the default is local `data/`.
Generated results and the manuscript are excluded from the software archive.

## Resume and individual stages

Rerun the same command to reuse completed results with matching processing
parameters. An unfinished fresh repeat reacquires its input. SRA, FASTQ, and BAM
intermediates are removed after their dependent results and input hashes are saved.
Use `-KeepFastq` or `-KeepBam` on `run_sample.ps1` or `run_repeat.ps1` to retain them.

```powershell
pwsh -NoProfile -File scripts/prepare_reference.ps1 -HostData D:/RNAseq-data
pwsh -NoProfile -File scripts/process_samples.ps1 -HostData D:/RNAseq-data
pwsh -NoProfile -File scripts/run_repeat.ps1 -Layout PAIRED -HostData D:/RNAseq-data -RunStrandBenchmark
pwsh -NoProfile -File scripts/run_repeat.ps1 -Layout SINGLE -HostData D:/RNAseq-data
python scripts/build_report.py --data-root D:/RNAseq-data
```

The complete report requires both repeats and the strand benchmark. Change result,
work, and log namespaces in `config.yaml` when changing scientific parameters.

## Inputs and QC

`config.yaml` fixes reference hashes, container digests, resources, and processing
parameters: STAR `Basic`/`intronMotif`, featureCounts `-s 2`, and RegTools PE=`RF`,
SE=`XS`. `requirements.lock` fixes Python wheels for CPython 3.13 on Windows x86-64.

`samples.tsv` supplies the processing manifest; `source_metadata.tsv` records
GEO/SRA metadata retrieved on 2026-09-26 from studies GSE190899 and GSE184456.
Processing follows the deposited SRA layout. Nominal lengths are PE100/SE75;
FastQC-observed lengths are 100 nt for PE and 35–76 nt for SE. Length fields
describe inputs and do not constrain FASTQ processing.

Review FastQC/MultiQC for each execution. The published protocol retains reads
without automatic trimming; record QC warnings and review decisions with that
run. The marker panel is descriptive, with one library per condition.

Code and configuration: MIT. Documentation: CC BY 4.0. See `LICENSE`.
Citation metadata are in `CITATION.cff`.
