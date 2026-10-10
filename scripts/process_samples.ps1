param(
  [string]$ConfigPath = '',
  [string]$HostData = '',
  [int]$Threads = 0,
  [double]$MinFreeGb = 0,
  [switch]$KeepFastq,
  [switch]$KeepBam
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'common.ps1')
$BranchDir = Split-Path -Parent $PSScriptRoot
$ConfigPath = Get-CanonicalConfigPath $ConfigPath
$Config = Get-BranchConfig $ConfigPath
$HostData = Get-DataRoot $Config $HostData
if ($Threads -le 0) { $Threads = [int]$Config.resources.threads }
if ($MinFreeGb -le 0) { $MinFreeGb = [double]$Config.resources.minimum_free_gb }
$SamplesPath = Join-Path $BranchDir 'samples.tsv'
$ResultRoot = Get-DataPath $HostData $Config.paths.data_results
$LogRoot = Get-DataPath $HostData $Config.paths.data_logs
$ResultNamespace = Split-Path -Leaf $ResultRoot
$ReferenceReady = Join-Path $HostData 'reference/BRANCH_A_REFERENCE_READY.json'
$RseqcBed = Join-Path $HostData "reference/$($Config.library.rseqc_bed)"
$Gtf = Join-Path $HostData "reference/$($Config.reference.gtf)"
$Runner = Join-Path $PSScriptRoot 'run_sample.ps1'
$StatusPath = Join-Path $LogRoot 'run_status.json'

New-Item -ItemType Directory -Force -Path $ResultRoot, $LogRoot | Out-Null

function Write-RunStatus([string]$Phase, [hashtable]$Extra = @{}) {
  $payload = @{ phase = $Phase; time = (Get-Date).ToString('o') }
  foreach ($key in $Extra.Keys) { $payload[$key] = $Extra[$key] }
  $payload | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $StatusPath -Encoding utf8
}

if (-not (Test-Path -LiteralPath $ReferenceReady)) {
  throw "The methods reference is not ready. Run prepare_reference.ps1 first: $ReferenceReady"
}
if (-not (Test-Path -LiteralPath $RseqcBed) -or (Get-Item -LiteralPath $RseqcBed).Length -lt 1MB) {
  Write-RunStatus 'build_rseqc_bed'
  $buildingBed = "$RseqcBed.building"
  Remove-Item -LiteralPath $buildingBed -Force -ErrorAction SilentlyContinue
  & python (Join-Path $PSScriptRoot 'gtf_to_rseqc_bed.py') --gtf $Gtf --bed $buildingBed
  if ($LASTEXITCODE -ne 0 -or -not (Test-Path -LiteralPath $buildingBed) -or (Get-Item -LiteralPath $buildingBed).Length -lt 1MB) { throw 'Failed to create RSeQC BED' }
  Move-Item -LiteralPath $buildingBed -Destination $RseqcBed -Force
}

$samples = Import-Csv -LiteralPath $SamplesPath -Delimiter "`t"
$completed = New-Object System.Collections.Generic.List[string]
foreach ($sample in $samples) {
  Write-RunStatus 'sample' @{ sample_id = $sample.sample_id; completed = @($completed) }
  $args = @{
    ConfigPath = $ConfigPath
    SampleId = $sample.sample_id
    Srr = $sample.SRR
    Layout = $sample.layout
    HostData = $HostData
    Threads = $Threads
    MinFreeGb = $MinFreeGb
  }
  if ($KeepFastq) { $args.KeepFastq = $true }
  if ($KeepBam) { $args.KeepBam = $true }
  & $Runner @args
  if ($LASTEXITCODE -ne 0) { throw "Sample runner failed for $($sample.sample_id)" }
  $completed.Add($sample.sample_id)
}

Write-RunStatus 'multiqc' @{ completed = @($completed) }
$multiqcDir = Join-Path $ResultRoot 'multiqc'
New-Item -ItemType Directory -Force -Path $multiqcDir | Out-Null
$multiqcConfig = Join-Path $ResultRoot 'multiqc_config.yaml'
Copy-Item -LiteralPath (Join-Path $BranchDir 'multiqc_config.yaml') -Destination $multiqcConfig -Force
$multiqcLog = Join-Path $LogRoot 'multiqc.log'
$multiqcImage = $Config.containers.multiqc
$multiqcMemory = [double]$Config.resources.stage_memory_gb.multiqc
$qcInputs = @(@('fastqc', 'counts', 'align_qc', 'rseqc') | ForEach-Object { "/data/results/$ResultNamespace/$_" })
# Use only the completed second-pass STAR logs, excluding pass-1 scratch logs.
$qcInputs += @($samples | ForEach-Object { "/data/results/$ResultNamespace/star/$($_.sample_id)/$($_.sample_id).Log.final.out" })
$multiqcInputs = $qcInputs -join ' '
$multiqcCmd = "docker run --rm --memory=${multiqcMemory}g --cpus=2 -v `"${HostData}:/data`" -w /data $multiqcImage multiqc $multiqcInputs --config /data/results/$ResultNamespace/multiqc_config.yaml --filename multiqc_report.html --force --outdir /data/results/$ResultNamespace/multiqc"
cmd.exe /d /s /c "$multiqcCmd 1> `"$multiqcLog`" 2>&1"
if ($LASTEXITCODE -ne 0) { throw "MultiQC failed; see $multiqcLog" }

Write-RunStatus 'done' @{ completed = @($completed); sample_count = $completed.Count }
Write-Host "Methods workflow complete: $($completed.Count) samples"
