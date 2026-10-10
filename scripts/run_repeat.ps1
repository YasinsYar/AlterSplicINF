param(
  [string]$ConfigPath = '',
  [string]$HostData = '',
  [ValidateSet('PAIRED', 'SINGLE')][string]$Layout,
  [int]$Threads = 0,
  [switch]$RunStrandBenchmark,
  [switch]$KeepBam,
  [switch]$KeepFastq
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'common.ps1')
$ConfigPath = Get-CanonicalConfigPath $ConfigPath
$Config = Get-BranchConfig $ConfigPath
$BranchDir = Get-BranchRoot
$HostData = Get-DataRoot $Config $HostData
if ($Threads -le 0) { $Threads = [int]$Config.resources.threads }
$RepeatConfig = if ($Layout -eq 'PAIRED') { $Config.reprocheck } else { $Config.reprocheck_se }
$MonitorInterval = [int]$Config.resources.monitor_interval_seconds
$Sample = $RepeatConfig.sample_id
$Srr = $RepeatConfig.srr
$Namespace = $RepeatConfig.result_namespace
$ResultRoot = Join-Path $HostData "results/$Namespace"
$WorkRoot = Join-Path $HostData "work/$Namespace"
$LogRoot = Join-Path (Get-DataPath $HostData $Config.paths.data_logs) $RepeatConfig.log_namespace
$StarDir = Join-Path $ResultRoot 'star'
$FastqcDir = Join-Path $ResultRoot 'fastqc'
$AlignQcDir = Join-Path $ResultRoot 'align_qc'
$CountDir = Join-Path $ResultRoot 'counts'
$JunctionDir = Join-Path $ResultRoot 'junctions'
$RseqcDir = Join-Path $ResultRoot 'rseqc'
$MetricDir = Join-Path $ResultRoot 'metrics'
$StatusPath = Join-Path $LogRoot 'status.json'
$CompletePath = Join-Path $ResultRoot 'REPROCHECK_COMPLETE.json'
$MetricPath = Join-Path $MetricDir 'fresh_end_to_end.tsv'
$CommandPath = Join-Path $LogRoot 'commands.jsonl'
$InputManifestPath = Join-Path $ResultRoot 'input_manifest.json'
$StageMemory = $Config.resources.stage_memory_gb
$Repeat = "${Sample}_fresh_repeat"
$BenchmarkDir = Join-Path $ResultRoot 'strand_benchmark'
$BenchmarkMetricPath = Join-Path $MetricDir 'strand_benchmark.tsv'
$PrimaryInputPath = Join-Path (Get-DataPath $HostData $Config.paths.data_results) "input_checksums/$Sample.json"
if (-not (Test-Path -LiteralPath $PrimaryInputPath)) { throw 'Primary input manifest is missing; run process_samples.ps1 first.' }
$RegtoolsStrand = Get-RegtoolsStrand $Config $Layout
$JunctionSignature = Get-JunctionPolicySignature $Config $Layout
$ProcessingSignature = Get-ProcessingSignature $Config $Srr $Layout

foreach ($dir in @($ResultRoot, $WorkRoot, $LogRoot, $StarDir, $FastqcDir, $AlignQcDir, $CountDir, $JunctionDir, $RseqcDir, $MetricDir)) {
  New-Item -ItemType Directory -Force -Path $dir | Out-Null
}

function Get-FreeGb { Get-FreeSpaceGb $HostData }

function Get-PathBytes([string[]]$Paths) {
  [int64]$total = 0
  foreach ($path in $Paths) {
    if (Test-Path -LiteralPath $path -PathType Leaf) { $total += (Get-Item -LiteralPath $path).Length }
    elseif (Test-Path -LiteralPath $path -PathType Container) {
      $sum = (Get-ChildItem -LiteralPath $path -File -Recurse -ErrorAction SilentlyContinue | Measure-Object Length -Sum).Sum
      if ($null -ne $sum) { $total += $sum }
    }
  }
  return $total
}

function Write-Status([string]$Phase, [hashtable]$Extra = @{}) {
  $payload = @{ phase = $Phase; sample_id = $Sample; time = (Get-Date).ToString('o'); free_gb = Get-FreeGb }
  foreach ($key in $Extra.Keys) { $payload[$key] = $Extra[$key] }
  $payload | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $StatusPath -Encoding utf8
}

function Invoke-MeasuredStep {
  param([string]$Name, [string]$Command, [string[]]$OutputPaths = @())
  $log = Join-Path $LogRoot "$Name.log"
  $started = Get-Date
  $freeBefore = Get-FreeGb
  Write-Status $Name @{ command = $Command }
  @{
    stage = $Name; started_at = $started.ToString('o'); command = $Command;
    config_sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $ConfigPath).Hash.ToLowerInvariant()
  } | ConvertTo-Json -Compress | Add-Content -LiteralPath $CommandPath -Encoding utf8
  $monitor = Start-Job -ArgumentList $MonitorInterval, $HostData -ScriptBlock {
    param($Seconds, $DataPath)
    while ($true) {
      $root = [IO.Path]::GetPathRoot([IO.Path]::GetFullPath($DataPath))
      $driveName = $root.TrimEnd([IO.Path]::DirectorySeparatorChar, [IO.Path]::AltDirectorySeparatorChar).TrimEnd(':')
      $drive = Get-PSDrive -Name $driveName
      [math]::Round($drive.Free / 1GB, 3)
      Start-Sleep -Seconds $Seconds
    }
  }
  cmd.exe /d /s /c "$Command 1> `"$log`" 2>&1"
  $exitCode = $LASTEXITCODE
  Stop-Job -Job $monitor -ErrorAction SilentlyContinue
  $freeSamples = @($freeBefore) + @(Receive-Job -Job $monitor -ErrorAction SilentlyContinue)
  Remove-Job -Job $monitor -Force -ErrorAction SilentlyContinue
  $finished = Get-Date
  $freeAfter = Get-FreeGb
  $minimumFree = ($freeSamples | Measure-Object -Minimum).Minimum
  if ($null -eq $minimumFree) { $minimumFree = [math]::Min($freeBefore, $freeAfter) }
  if (-not (Test-Path -LiteralPath $MetricPath)) {
    "stage`tstarted_at`tfinished_at`telapsed_seconds`tfree_before_gb`tminimum_free_gb`tfree_after_gb`ttransient_disk_delta_gb`toutput_bytes`texit_code" | Set-Content -LiteralPath $MetricPath -Encoding utf8
  }
  $culture = [Globalization.CultureInfo]::InvariantCulture
  @(
    $Name, $started.ToString('o'), $finished.ToString('o'), (($finished-$started).TotalSeconds).ToString('0.0',$culture),
    ([double]$freeBefore).ToString('0.000',$culture), ([double]$minimumFree).ToString('0.000',$culture),
    ([double]$freeAfter).ToString('0.000',$culture), ([math]::Max(0,$freeBefore-$minimumFree)).ToString('0.000',$culture),
    (Get-PathBytes $OutputPaths).ToString($culture), $exitCode.ToString($culture)
  ) -join "`t" | Add-Content -LiteralPath $MetricPath -Encoding utf8
  if ($exitCode -ne 0) {
    Write-Status 'failed' @{ step = $Name; exit_code = $exitCode; log = $log }
    throw "Fresh repeat step $Name failed; see $log"
  }
}

if (Test-Path -LiteralPath $CompletePath) {
  $prior = Get-Content -LiteralPath $CompletePath -Raw | ConvertFrom-Json
  if (-not $prior.ready -or $prior.junction_policy_signature -ne $JunctionSignature) {
    throw 'The existing repeat uses another policy. Choose a new result namespace.'
  }
  if ($prior.processing_signature -and $prior.processing_signature -ne $ProcessingSignature) { throw 'Existing repeat processing policy differs; choose a new namespace.' }
  if ($RunStrandBenchmark -and -not (Test-Path -LiteralPath $BenchmarkMetricPath)) {
    throw 'The completed repeat has no strand experiment. Choose a new paired repeat namespace and use -RunStrandBenchmark.'
  }
  foreach ($required in @($InputManifestPath, (Join-Path $CountDir "$Repeat.featureCounts.txt"),
    (Join-Path $CountDir "$Repeat.featureCounts.txt.summary"), (Join-Path $JunctionDir "$Repeat.junc"),
    (Join-Path $RseqcDir "$Repeat.infer_experiment.txt"), (Join-Path $StarDir "$Repeat.Log.final.out"), $MetricPath)) {
    if (-not (Test-Path -LiteralPath $required)) { throw "Completed repeat artifact is missing: $required" }
  }
  Write-Host "$Layout repeat already complete."
  exit 0
}
if ((Get-FreeGb) -lt [double]$Config.resources.minimum_free_gb) { throw 'Insufficient free disk for fresh repeat' }

if (Test-Path -LiteralPath $MetricPath) { Remove-Item -LiteralPath $MetricPath -Force }
if (Test-Path -LiteralPath $CommandPath) { Remove-Item -LiteralPath $CommandPath -Force }
$sraDir = Join-Path $WorkRoot "sra/$Sample"
$fastqDir = Join-Path $WorkRoot "fastq/$Sample"
Remove-WorkDirectory $sraDir $WorkRoot
Remove-WorkDirectory $fastqDir $WorkRoot
New-Item -ItemType Directory -Force -Path $sraDir, $fastqDir | Out-Null

$sraImage = $Config.containers.sra_tools
$pigzImage = $Config.containers.pigz
$starImage = $Config.containers.star
$fastqcImage = $Config.containers.fastqc
$samtoolsImage = $Config.containers.samtools
$subreadImage = $Config.containers.subread
$regtoolsImage = $Config.containers.regtools
$rseqcImage = $Config.containers.rseqc
$freshR1 = Join-Path $fastqDir "${Sample}_R1.fastq.gz"
$freshR2 = Join-Path $fastqDir "${Sample}_R2.fastq.gz"

$prefetch = "docker run --rm --name altersplcinf-a-se-prefetch --memory=$($StageMemory.prefetch)g --cpus=$Threads -v `"${HostData}:/data`" -w /data $sraImage prefetch $Srr -O /data/work/$Namespace/sra/$Sample"
Invoke-MeasuredStep 'prefetch' $prefetch @($sraDir)
$sraFile = Get-ChildItem -LiteralPath $sraDir -Filter "$Srr.sra" -Recurse -File | Select-Object -First 1
if ($null -eq $sraFile) { throw "Fresh prefetch did not produce $Srr.sra" }
$sraSha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $sraFile.FullName).Hash.ToLowerInvariant()
$relativeSra = $sraFile.FullName.Substring($HostData.Replace('/', '\').TrimEnd('\').Length).Replace('\','/')
$fasterq = "docker run --rm --name altersplcinf-a-se-fasterq --memory=$($StageMemory.fasterq)g --cpus=$Threads -v `"${HostData}:/data`" -w /data $sraImage fasterq-dump /data$relativeSra --split-files --threads $Threads -O /data/work/$Namespace/fastq/$Sample"
Invoke-MeasuredStep 'fasterq' $fasterq @($fastqDir)
$compressScript = "set -euo pipefail; d=/data/work/$Namespace/fastq/$Sample; if [ -f `"`$d/${Srr}_1.fastq`" ]; then pigz -p $Threads `"`$d/${Srr}_1.fastq`"; mv `"`$d/${Srr}_1.fastq.gz`" `"`$d/${Sample}_R1.fastq.gz`"; elif [ -f `"`$d/${Srr}.fastq`" ]; then pigz -p $Threads `"`$d/${Srr}.fastq`"; mv `"`$d/${Srr}.fastq.gz`" `"`$d/${Sample}_R1.fastq.gz`"; else echo FASTQ_missing >&2; exit 1; fi"
if ($Layout -eq 'PAIRED') {
  $compressScript = "set -euo pipefail; d=/data/work/$Namespace/fastq/$Sample; pigz -p $Threads `$d/${Srr}_1.fastq; pigz -p $Threads `$d/${Srr}_2.fastq; mv `$d/${Srr}_1.fastq.gz `$d/${Sample}_R1.fastq.gz; mv `$d/${Srr}_2.fastq.gz `$d/${Sample}_R2.fastq.gz"
}
$compress = "docker run --rm --name altersplcinf-a-se-compress --memory=$($StageMemory.compression)g --cpus=$Threads -v `"${HostData}:/data`" -w /data $pigzImage bash -lc `"$compressScript`""
Invoke-MeasuredStep 'compress_fastq' $compress @($freshR1)
$hashPaths = @($freshR1)
if ($Layout -eq 'PAIRED') { $hashPaths += $freshR2 }
$fastqJson = & python (Join-Path $PSScriptRoot 'hash_fastq_content.py') @hashPaths
if ($LASTEXITCODE -ne 0) { throw 'Repeat FASTQ hashing failed.' }
$records = @((($fastqJson | Out-String | ConvertFrom-Json).PSObject.Properties.Value))
$inputManifest = [ordered]@{ sample_id = $Sample; srr = $Srr; layout = $Layout;
  sra_size_bytes = $sraFile.Length; sra_sha256 = $sraSha256; fastq = $records[0] }
if ($Layout -eq 'PAIRED') { $inputManifest.fastq_r2 = $records[1] }
$inputManifest | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $InputManifestPath -Encoding utf8
$primary = Get-Content -LiteralPath $PrimaryInputPath -Raw | ConvertFrom-Json
foreach ($mate in $(if ($Layout -eq 'PAIRED') { @('fastq', 'fastq_r2') } else { @('fastq') })) {
  if ($primary.$mate.decompressed_content_sha256 -ne $inputManifest[$mate].decompressed_content_sha256) {
    throw "Repeat FASTQ content differs from primary input ($mate)."
  }
}
Remove-WorkDirectory $sraDir $WorkRoot
$readFiles = "/data/work/$Namespace/fastq/$Sample/${Sample}_R1.fastq.gz"
if ($Layout -eq 'PAIRED') { $readFiles += " /data/work/$Namespace/fastq/$Sample/${Sample}_R2.fastq.gz" }

$fastqc = "docker run --rm --name altersplcinf-a-se-fastqc --memory=$($StageMemory.fastqc)g --cpus=2 -v `"${HostData}:/data`" -w /data $fastqcImage fastqc -t 2 -o /data/results/$Namespace/fastqc $readFiles"
Invoke-MeasuredStep 'fastqc' $fastqc @($FastqcDir)
$bam = "/data/results/$Namespace/star/$Repeat.Aligned.sortedByCoord.out.bam"
$star = "docker run --rm --name altersplcinf-a-se-star --memory=$($StageMemory.star)g --cpus=$Threads -v `"${HostData}:/data`" -w /data $starImage STAR --genomeDir /data/reference/$($Config.reference.star_index) --readFilesIn $readFiles --readFilesCommand zcat --runThreadN $Threads --twopassMode $($Config.library.star_twopass_mode) --outSAMtype BAM SortedByCoordinate --outSAMstrandField $($Config.library.star_out_sam_strand_field) --outFileNamePrefix /data/results/$Namespace/star/$Repeat."
Invoke-MeasuredStep 'star_2pass' $star @($StarDir)
if (-not (Test-CompleteBam (Join-Path $StarDir "$Repeat.Aligned.sortedByCoord.out.bam"))) { throw 'STAR BAM is incomplete.' }
$samtoolsIndex = "docker run --rm --name altersplcinf-a-se-samtools-index --memory=$($StageMemory.samtools)g --cpus=2 -v `"${HostData}:/data`" -w /data $samtoolsImage samtools index -@ 2 $bam"
Invoke-MeasuredStep 'samtools_index' $samtoolsIndex @((Join-Path $StarDir "$Repeat.Aligned.sortedByCoord.out.bam.bai"))
$flagstat = "docker run --rm --name altersplcinf-a-se-samtools-flagstat --memory=$($StageMemory.samtools)g --cpus=2 -v `"${HostData}:/data`" -w /data $samtoolsImage bash -lc `"samtools flagstat -@ 2 '$bam' > '/data/results/$Namespace/align_qc/$Repeat.flagstat.txt'`""
Invoke-MeasuredStep 'samtools_flagstat' $flagstat @((Join-Path $AlignQcDir "$Repeat.flagstat.txt"))
$pairOptions = if ($Layout -eq 'PAIRED') { '-p --countReadPairs' } else { '' }
$count = "docker run --rm --name altersplcinf-a-se-featurecounts --memory=$($StageMemory.featurecounts)g --cpus=$Threads -v `"${HostData}:/data`" -w /data $subreadImage featureCounts -T $Threads -a /data/reference/$($Config.reference.gtf) -o /data/results/$Namespace/counts/$Repeat.featureCounts.txt -s $($Config.library.featurecounts_s) $pairOptions $bam"
Invoke-MeasuredStep 'featurecounts_s2' $count @((Join-Path $CountDir "$Repeat.featureCounts.txt"),(Join-Path $CountDir "$Repeat.featureCounts.txt.summary"))
$junction = "docker run --rm --name altersplcinf-a-se-regtools --memory=$($StageMemory.regtools)g --cpus=2 -v `"${HostData}:/data`" -w /data $regtoolsImage regtools junctions extract -a $($Config.library.regtools_anchor) -m $($Config.library.regtools_min_intron) -M $($Config.library.regtools_max_intron) -s $RegtoolsStrand $bam -o /data/results/$Namespace/junctions/$Repeat.junc"
Invoke-MeasuredStep "regtools_$($RegtoolsStrand.ToLowerInvariant())" $junction @((Join-Path $JunctionDir "$Repeat.junc"))
$checkArgs = @((Join-Path $PSScriptRoot 'junction_strand_check.py'), '--junctions', (Join-Path $JunctionDir "$Repeat.junc"), '--report', (Join-Path $JunctionDir "$Repeat.strand_check.json"))
if ($Layout -eq 'SINGLE') { $checkArgs += '--require-both-strands' }
& python @checkArgs
if ($LASTEXITCODE -ne 0) { throw 'Repeat junction strand check failed; BAM is retained.' }
$rseqc = "docker run --rm --name altersplcinf-a-se-rseqc --memory=$($StageMemory.rseqc)g --cpus=2 -v `"${HostData}:/data`" -w /data $rseqcImage bash -lc `"infer_experiment.py -r /data/reference/$($Config.library.rseqc_bed) -i '$bam' > '/data/results/$Namespace/rseqc/$Repeat.infer_experiment.txt'`""
Invoke-MeasuredStep 'rseqc' $rseqc @((Join-Path $RseqcDir "$Repeat.infer_experiment.txt"))

if ($RunStrandBenchmark) {
  if ($Layout -ne 'PAIRED') { throw 'The controlled strand experiment is defined for the paired representative.' }
  New-Item -ItemType Directory -Force -Path $BenchmarkDir | Out-Null
  Copy-Item -LiteralPath (Join-Path $CountDir "$Repeat.featureCounts.txt.summary") -Destination (Join-Path $BenchmarkDir "$Repeat.s2.featureCounts.txt.summary") -Force
  $productionMetricPath = $MetricPath
  $MetricPath = $BenchmarkMetricPath
  Remove-Item -LiteralPath $MetricPath -Force -ErrorAction SilentlyContinue
  foreach ($mode in @(0, 1)) {
    $output = "/data/results/$Namespace/strand_benchmark/$Repeat.s${mode}.featureCounts.txt"
    $command = "docker run --rm --name altersplcinf-a-strand-s$mode --memory=$($StageMemory.featurecounts)g --cpus=$Threads -v `"${HostData}:/data`" -w /data $subreadImage featureCounts -T $Threads -a /data/reference/$($Config.reference.gtf) -o $output -s $mode -p --countReadPairs $bam"
    Invoke-MeasuredStep "featurecounts_s$mode" $command @((Join-Path $BenchmarkDir "$Repeat.s${mode}.featureCounts.txt.summary"))
    Remove-Item -LiteralPath (Join-Path $BenchmarkDir "$Repeat.s${mode}.featureCounts.txt") -Force
  }
  $MetricPath = $productionMetricPath
}

$bamHost = Join-Path $StarDir "$Repeat.Aligned.sortedByCoord.out.bam"
if (-not $KeepBam) { Remove-Item -LiteralPath $bamHost, "$bamHost.bai" -Force -ErrorAction SilentlyContinue }
if (-not $KeepFastq) { Remove-WorkDirectory $fastqDir $WorkRoot }
$complete = @{
  ready = $true; repeat_type = 'fresh_from_sra'; sample_id = $Sample; srr = $Srr; layout = $Layout;
  completed = (Get-Date).ToString('o'); config_sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $ConfigPath).Hash.ToLowerInvariant();
  metrics = $MetricPath; input_manifest = $InputManifestPath; free_gb = Get-FreeGb;
  junction_policy_signature = $JunctionSignature; regtools_strand = $RegtoolsStrand; processing_signature = $ProcessingSignature
}
$complete | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $CompletePath -Encoding utf8
Write-Status 'done' @{ complete_manifest = $CompletePath }
