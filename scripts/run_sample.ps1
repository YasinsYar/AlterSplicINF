param(
  [Parameter(Mandatory = $true)][string]$SampleId,
  [Parameter(Mandatory = $true)][string]$Srr,
  [Parameter(Mandatory = $true)][ValidateSet('PAIRED', 'SINGLE')][string]$Layout,
  [string]$ConfigPath = '',
  [string]$HostData = '',
  [int]$Threads = 0,
  [double]$MinFreeGb = 0,
  [switch]$KeepFastq,
  [switch]$KeepBam
)

$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'
. (Join-Path $PSScriptRoot 'common.ps1')
$ConfigPath = Get-CanonicalConfigPath $ConfigPath
$Config = Get-BranchConfig $ConfigPath
$HostData = Get-DataRoot $Config $HostData
if ($Threads -le 0) { $Threads = [int]$Config.resources.threads }
if ($MinFreeGb -le 0) { $MinFreeGb = [double]$Config.resources.minimum_free_gb }

$Images = @{
  sra = $Config.containers.sra_tools
  pigz = $Config.containers.pigz
  fastqc = $Config.containers.fastqc
  star = $Config.containers.star
  samtools = $Config.containers.samtools
  subread = $Config.containers.subread
  regtools = $Config.containers.regtools
  rseqc = $Config.containers.rseqc
}

$ResultRoot = Get-DataPath $HostData $Config.paths.data_results
$WorkRoot = Get-DataPath $HostData $Config.paths.data_work
$ResultNamespace = Split-Path -Leaf $ResultRoot
$WorkNamespace = Split-Path -Leaf $WorkRoot
$LogRoot = Join-Path (Get-DataPath $HostData $Config.paths.data_logs) "pipeline/samples/$SampleId"
$StarDir = Join-Path $ResultRoot "star/$SampleId"
$FastqcDir = Join-Path $ResultRoot 'fastqc'
$AlignQcDir = Join-Path $ResultRoot 'align_qc'
$CountDir = Join-Path $ResultRoot 'counts'
$JunctionDir = Join-Path $ResultRoot 'junctions'
$RseqcDir = Join-Path $ResultRoot 'rseqc'
$MetricDir = Join-Path $ResultRoot 'metrics/pipeline'
$MetricPath = Join-Path $MetricDir "$SampleId.steps.tsv"
$CommandManifestPath = Join-Path $LogRoot 'commands.jsonl'
$StatusPath = Join-Path $LogRoot 'status.json'
$Fasta = Join-Path $HostData ("reference/" + $Config.reference.fasta)
$Gtf = Join-Path $HostData ("reference/" + $Config.reference.gtf)
$Index = Join-Path $HostData ("reference/" + $Config.reference.star_index)
$RseqcBed = Join-Path $HostData ("reference/" + $Config.library.rseqc_bed)
$StageMemory = $Config.resources.stage_memory_gb
$RegtoolsStrand = Get-RegtoolsStrand $Config $Layout
$JunctionSignature = Get-JunctionPolicySignature $Config $Layout
$ProcessingSignature = Get-ProcessingSignature $Config $Srr $Layout
$JunctionPolicyPath = Join-Path $JunctionDir "$SampleId.junc.policy.json"
$junctionPolicyMatches = $false
if (Test-Path -LiteralPath $JunctionPolicyPath) {
  $savedPolicy = Get-Content -LiteralPath $JunctionPolicyPath -Raw | ConvertFrom-Json
  $junctionPolicyMatches = $savedPolicy.signature -eq $JunctionSignature -and $savedPolicy.srr -eq $Srr
}

foreach ($dir in @($ResultRoot, $WorkRoot, $LogRoot, $StarDir, $FastqcDir, $AlignQcDir, $CountDir, $JunctionDir, $RseqcDir, $MetricDir)) {
  New-Item -ItemType Directory -Force -Path $dir | Out-Null
}

foreach ($required in @($Fasta, $Gtf, (Join-Path $Index 'Genome'), $RseqcBed)) {
  if (-not (Test-Path -LiteralPath $required)) { throw "Missing methods-workflow prerequisite: $required" }
}

function Get-FreeGb {
  Get-FreeSpaceGb $HostData
}

function Get-PathBytes([string[]]$Paths) {
  [int64]$total = 0
  foreach ($path in $Paths) {
    if (Test-Path -LiteralPath $path -PathType Leaf) {
      $total += (Get-Item -LiteralPath $path).Length
    } elseif (Test-Path -LiteralPath $path -PathType Container) {
      $sum = (Get-ChildItem -LiteralPath $path -File -Recurse -ErrorAction SilentlyContinue | Measure-Object Length -Sum).Sum
      if ($null -ne $sum) { $total += $sum }
    }
  }
  return $total
}

function Convert-MemoryToMb([string]$Value) {
  if ($Value -notmatch '^\s*([0-9.]+)\s*([KMG]i?B)') { return 0.0 }
  $number = [double]::Parse($Matches[1], [Globalization.CultureInfo]::InvariantCulture)
  switch -Regex ($Matches[2]) {
    '^Ki?B$' { return $number / 1024 }
    '^Gi?B$' { return $number * 1024 }
    default { return $number }
  }
}

function Write-Status([string]$Phase, [hashtable]$Extra = @{}) {
  $payload = @{
    sample_id = $SampleId
    srr = $Srr
    phase = $Phase
    processing_signature = $ProcessingSignature
    time = (Get-Date).ToString('o')
    free_gb = Get-FreeGb
  }
  foreach ($key in $Extra.Keys) { $payload[$key] = $Extra[$key] }
  $payload | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $StatusPath -Encoding utf8
}

function Add-Metric {
  param(
    [string]$Stage,
    [datetime]$Started,
    [datetime]$Finished,
    [int]$CpuLimit,
    [double]$MemoryLimitGb,
    [double]$PeakMemoryMb,
    [double]$PeakCpuPercent,
    [int64]$OutputBytes,
    [int]$ExitCode,
    [string]$State
  )
  if (-not (Test-Path -LiteralPath $MetricPath)) {
    "sample_id`tstage`tstarted_at`tfinished_at`telapsed_seconds`tcpu_limit`tmemory_limit_gb`tpeak_memory_mb`tpeak_cpu_percent`toutput_bytes`texit_code`tstate" |
      Set-Content -LiteralPath $MetricPath -Encoding utf8
  }
  $elapsed = [math]::Round(($Finished - $Started).TotalSeconds, 1)
  $culture = [Globalization.CultureInfo]::InvariantCulture
  $line = @(
    $SampleId, $Stage, $Started.ToString('o'), $Finished.ToString('o'), $elapsed.ToString('0.0', $culture),
    $CpuLimit.ToString($culture), $MemoryLimitGb.ToString('0.###', $culture), ([math]::Round($PeakMemoryMb, 1)).ToString('0.0', $culture),
    ([math]::Round($PeakCpuPercent, 1)).ToString('0.0', $culture), $OutputBytes.ToString($culture), $ExitCode.ToString($culture), $State
  ) -join "`t"
  Add-Content -LiteralPath $MetricPath -Value $line -Encoding utf8
}

function Quote-CommandArgument([string]$Value) {
  if ($Value -eq '') { return '""' }
  if ($Value -match '[\s"&|<>^]') { return '"' + ($Value -replace '"', '\"') + '"' }
  return $Value
}

function Invoke-MeasuredDocker {
  param(
    [Parameter(Mandatory = $true)][string]$Stage,
    [Parameter(Mandatory = $true)][string]$Image,
    [Parameter(Mandatory = $true)][string[]]$Arguments,
    [Parameter(Mandatory = $true)][string]$LogPath,
    [Parameter(Mandatory = $true)][double]$MemoryGb,
    [int]$Cpu = $Threads,
    [string[]]$OutputPaths = @()
  )
  $safeStage = $Stage.ToLowerInvariant() -replace '[^a-z0-9]+', '-'
  $safeSample = $SampleId.ToLowerInvariant() -replace '[^a-z0-9]+', '-'
  $container = "altersplcinf-a-$safeSample-$safeStage-$PID"
  $dockerArgs = @('run', '--rm', '--name', $container, "--memory=${MemoryGb}g", "--cpus=$Cpu", '-v', "${HostData}:/data", '-w', '/data', $Image) + $Arguments
  $command = 'docker ' + (($dockerArgs | ForEach-Object { Quote-CommandArgument $_ }) -join ' ')
  $started = Get-Date
  Write-Status $Stage @{ command = $command }
  @{
    sample_id = $SampleId
    stage = $Stage
    time = $started.ToString('o')
    config_path = $ConfigPath
    config_sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $ConfigPath).Hash.ToLowerInvariant()
    command = $command
  } | ConvertTo-Json -Compress | Add-Content -LiteralPath $CommandManifestPath -Encoding utf8
  Write-Host "[$Stage] $command"

  $monitor = Start-Job -ScriptBlock {
    param($ContainerName)
    while ($true) {
      $value = & docker stats --no-stream --format '{{.MemUsage}}|{{.CPUPerc}}' $ContainerName 2>$null
      if ($value) { Write-Output $value }
      Start-Sleep -Seconds 2
    }
  } -ArgumentList $container

  cmd.exe /d /s /c "$command 1> `"$LogPath`" 2>&1"
  $exitCode = $LASTEXITCODE
  Stop-Job -Job $monitor -ErrorAction SilentlyContinue
  $samples = @(Receive-Job -Job $monitor -ErrorAction SilentlyContinue)
  Remove-Job -Job $monitor -Force -ErrorAction SilentlyContinue

  $peakMemory = 0.0
  $peakCpu = 0.0
  foreach ($sample in $samples) {
    $parts = "$sample".Split('|')
    if ($parts.Count -lt 2) { continue }
    $memory = Convert-MemoryToMb ($parts[0].Split('/')[0].Trim())
    $cpuText = $parts[1].Trim().TrimEnd('%')
    $cpuValue = 0.0
    [double]::TryParse($cpuText, [Globalization.NumberStyles]::Float, [Globalization.CultureInfo]::InvariantCulture, [ref]$cpuValue) | Out-Null
    if ($memory -gt $peakMemory) { $peakMemory = $memory }
    if ($cpuValue -gt $peakCpu) { $peakCpu = $cpuValue }
  }
  $finished = Get-Date
  $bytes = Get-PathBytes $OutputPaths
  Add-Metric -Stage $Stage -Started $started -Finished $finished -CpuLimit $Cpu -MemoryLimitGb $MemoryGb `
    -PeakMemoryMb $peakMemory -PeakCpuPercent $peakCpu -OutputBytes $bytes -ExitCode $exitCode -State 'executed'
  if ($exitCode -ne 0) {
    Write-Status 'failed' @{ step = $Stage; exit_code = $exitCode; log = $LogPath }
    throw "Docker step $Stage failed with exit $exitCode; see $LogPath"
  }
}

function Record-Reused([string]$Stage, [int]$Cpu, [double]$MemoryGb, [string[]]$OutputPaths) {
  $now = Get-Date
  Add-Metric -Stage $Stage -Started $now -Finished $now -CpuLimit $Cpu -MemoryLimitGb $MemoryGb `
    -PeakMemoryMb 0 -PeakCpuPercent 0 -OutputBytes (Get-PathBytes $OutputPaths) -ExitCode 0 -State 'reused'
  Write-Host "[$Stage] reuse existing output"
}

$finishedOutputs = @(
  (Join-Path $StarDir "${SampleId}.Log.final.out"),
  (Join-Path $CountDir "$SampleId.featureCounts.txt"),
  (Join-Path $CountDir "$SampleId.featureCounts.txt.summary"),
  (Join-Path $JunctionDir "$SampleId.junc"),
  (Join-Path $RseqcDir "$SampleId.infer_experiment.txt"),
  (Join-Path $FastqcDir "${SampleId}_R1_fastqc.zip")
)
if ($Layout -eq 'PAIRED') { $finishedOutputs += (Join-Path $FastqcDir "${SampleId}_R2_fastqc.zip") }
if ($SampleId -in @($Config.reprocheck.sample_id, $Config.reprocheck_se.sample_id)) { $finishedOutputs += (Join-Path $ResultRoot "input_checksums/$SampleId.json") }
$previousStatus = $null
if (Test-Path -LiteralPath $StatusPath) {
  $previousStatus = Get-Content -LiteralPath $StatusPath -Raw | ConvertFrom-Json
  if ($previousStatus.processing_signature -and $previousStatus.processing_signature -ne $ProcessingSignature) { throw 'Existing artifacts use another processing policy. Choose a new result/work/log namespace.' }
}
if ($null -ne $previousStatus -and $previousStatus.phase -eq 'done' -and $junctionPolicyMatches -and
    ($finishedOutputs | Where-Object { -not (Test-Path -LiteralPath $_) }).Count -eq 0) {
  Write-Host "[done] $SampleId already has a complete artifact set"
  exit 0
}
if ((Get-FreeGb) -lt $MinFreeGb) { throw "Disk gate failed: $(Get-FreeGb) GB free < $MinFreeGb GB" }

Write-Status 'starting'

$branchFastqDir = Join-Path $WorkRoot "fastq/$SampleId"
$fastqDir = $branchFastqDir
if (-not (Test-Path -LiteralPath (Join-Path $branchFastqDir "${SampleId}_R1.fastq.gz"))) {
  $sraDir = Join-Path $WorkRoot "sra/$SampleId"
  New-Item -ItemType Directory -Force -Path $sraDir, $branchFastqDir | Out-Null
  Invoke-MeasuredDocker -Stage 'prefetch' -Image $Images.sra -MemoryGb $StageMemory.prefetch -Cpu $Threads `
    -LogPath (Join-Path $LogRoot 'prefetch.log') `
    -Arguments @('prefetch', $Srr, '-O', "/data/work/$WorkNamespace/sra/$SampleId") `
    -OutputPaths @($sraDir)
  $sraFile = Get-ChildItem -LiteralPath $sraDir -Filter "$Srr.sra" -Recurse -File | Select-Object -First 1
  if ($null -eq $sraFile) { throw "prefetch completed but $Srr.sra was not found under $sraDir" }
  $relativeSra = $sraFile.FullName.Substring($HostData.Replace('/', '\').TrimEnd('\').Length).Replace('\', '/')
  Invoke-MeasuredDocker -Stage 'fasterq' -Image $Images.sra -MemoryGb $StageMemory.fasterq -Cpu $Threads `
    -LogPath (Join-Path $LogRoot 'fasterq.log') `
    -Arguments @('fasterq-dump', "/data$relativeSra", '--split-files', '--threads', "$Threads", '-O', "/data/work/$WorkNamespace/fastq/$SampleId") `
    -OutputPaths @($branchFastqDir)
  $pigzScript = @"
set -euo pipefail; d=/data/work/$WorkNamespace/fastq/$SampleId; if [ -f "`$d/${Srr}_1.fastq" ]; then pigz -p $Threads "`$d/${Srr}_1.fastq"; mv "`$d/${Srr}_1.fastq.gz" "`$d/${SampleId}_R1.fastq.gz"; if [ -f "`$d/${Srr}_2.fastq" ]; then pigz -p $Threads "`$d/${Srr}_2.fastq"; mv "`$d/${Srr}_2.fastq.gz" "`$d/${SampleId}_R2.fastq.gz"; fi; elif [ -f "`$d/${Srr}.fastq" ]; then pigz -p $Threads "`$d/${Srr}.fastq"; mv "`$d/${Srr}.fastq.gz" "`$d/${SampleId}_R1.fastq.gz"; else echo FASTQ_missing >&2; ls -la "`$d" >&2; exit 1; fi
"@
  Invoke-MeasuredDocker -Stage 'compress_fastq' -Image $Images.pigz -MemoryGb $StageMemory.compression -Cpu $Threads `
    -LogPath (Join-Path $LogRoot 'compress_fastq.log') -Arguments @('bash', '-lc', $pigzScript) `
    -OutputPaths @($branchFastqDir)
  Remove-WorkDirectory $sraDir $WorkRoot
} else {
  Record-Reused -Stage 'input_fastq' -Cpu 1 -MemoryGb 0 -OutputPaths @($branchFastqDir)
}

$hostR1 = Join-Path $fastqDir "${SampleId}_R1.fastq.gz"
$hostR2 = Join-Path $fastqDir "${SampleId}_R2.fastq.gz"
if (-not (Test-Path -LiteralPath $hostR1)) { throw "Missing R1 FASTQ: $hostR1" }
if ($Layout -eq 'PAIRED' -and -not (Test-Path -LiteralPath $hostR2)) { throw "Missing R2 FASTQ: $hostR2" }
$relativeFastq = $fastqDir.Substring($HostData.Replace('/', '\').TrimEnd('\').Length).Replace('\', '/')
$readPaths = @("/data$relativeFastq/${SampleId}_R1.fastq.gz")
if ($Layout -eq 'PAIRED') { $readPaths += "/data$relativeFastq/${SampleId}_R2.fastq.gz" }

$fastqcExpected = @((Join-Path $FastqcDir "${SampleId}_R1_fastqc.zip"))
if ($Layout -eq 'PAIRED') { $fastqcExpected += (Join-Path $FastqcDir "${SampleId}_R2_fastqc.zip") }
if (($fastqcExpected | Where-Object { -not (Test-Path -LiteralPath $_) }).Count -gt 0) {
  Invoke-MeasuredDocker -Stage 'fastqc' -Image $Images.fastqc -MemoryGb $StageMemory.fastqc -Cpu 2 `
    -LogPath (Join-Path $LogRoot 'fastqc.log') `
    -Arguments (@('fastqc', '-t', '2', '-o', "/data/results/$ResultNamespace/fastqc") + $readPaths) `
    -OutputPaths $fastqcExpected
} else { Record-Reused -Stage 'fastqc' -Cpu 2 -MemoryGb $StageMemory.fastqc -OutputPaths $fastqcExpected }

$bamHost = Join-Path $StarDir "${SampleId}.Aligned.sortedByCoord.out.bam"
$starLogHost = Join-Path $StarDir "${SampleId}.Log.final.out"
if (-not (Test-CompleteBam $bamHost) -or -not (Test-Path -LiteralPath $starLogHost)) {
  # An interrupted run can leave an empty BAM, an old final log, and scratch.
  # Clear only this sample's generated scratch before restarting STAR.
  $starFull = [IO.Path]::GetFullPath($StarDir).TrimEnd('\') + '\'
  foreach ($suffix in @('_STARtmp', '_STARpass1', '_STARgenome')) {
    $scratch = [IO.Path]::GetFullPath((Join-Path $StarDir "${SampleId}.$suffix"))
    if (-not $scratch.StartsWith($starFull, [StringComparison]::OrdinalIgnoreCase)) { throw 'Unsafe STAR scratch cleanup target.' }
    if (Test-Path -LiteralPath $scratch) { Remove-Item -LiteralPath $scratch -Recurse -Force }
  }
  $dependent = @("$bamHost.bai", (Join-Path $AlignQcDir "$SampleId.flagstat.txt"),
    (Join-Path $CountDir "$SampleId.featureCounts.txt"), (Join-Path $CountDir "$SampleId.featureCounts.txt.summary"),
    (Join-Path $JunctionDir "$SampleId.junc"), $JunctionPolicyPath,
    (Join-Path $RseqcDir "$SampleId.infer_experiment.txt"))
  foreach ($path in $dependent) { Remove-Item -LiteralPath $path -Force -ErrorAction SilentlyContinue }
  $junctionPolicyMatches = $false
  $starArgs = @('STAR', '--genomeDir', "/data/reference/$($Config.reference.star_index)", '--readFilesIn') + $readPaths + @(
    '--readFilesCommand', 'zcat', '--runThreadN', "$Threads", '--twopassMode', "$($Config.library.star_twopass_mode)",
    '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMstrandField', "$($Config.library.star_out_sam_strand_field)",
    '--outFileNamePrefix', "/data/results/$ResultNamespace/star/$SampleId/${SampleId}."
  )
  Invoke-MeasuredDocker -Stage 'star_2pass' -Image $Images.star -MemoryGb $StageMemory.star -Cpu $Threads `
    -LogPath (Join-Path $LogRoot 'star.log') -Arguments $starArgs -OutputPaths @($bamHost, $starLogHost)
} else { Record-Reused -Stage 'star_2pass' -Cpu $Threads -MemoryGb $StageMemory.star -OutputPaths @($bamHost, $starLogHost) }

$bamContainer = "/data/results/$ResultNamespace/star/$SampleId/${SampleId}.Aligned.sortedByCoord.out.bam"
$baiHost = "$bamHost.bai"
if (-not (Test-Path -LiteralPath $baiHost)) {
  Invoke-MeasuredDocker -Stage 'samtools_index' -Image $Images.samtools -MemoryGb $StageMemory.samtools -Cpu 2 `
    -LogPath (Join-Path $LogRoot 'samtools_index.log') `
    -Arguments @('samtools', 'index', '-@', '2', $bamContainer) -OutputPaths @($baiHost)
} else { Record-Reused -Stage 'samtools_index' -Cpu 2 -MemoryGb $StageMemory.samtools -OutputPaths @($baiHost) }

$flagstatHost = Join-Path $AlignQcDir "$SampleId.flagstat.txt"
if (-not (Test-Path -LiteralPath $flagstatHost)) {
  Invoke-MeasuredDocker -Stage 'samtools_flagstat' -Image $Images.samtools -MemoryGb $StageMemory.samtools -Cpu 2 `
    -LogPath (Join-Path $LogRoot 'samtools_flagstat.log') `
    -Arguments @('bash', '-lc', "samtools flagstat -@ 2 '$bamContainer' > '/data/results/$ResultNamespace/align_qc/$SampleId.flagstat.txt'") `
    -OutputPaths @($flagstatHost)
} else { Record-Reused -Stage 'samtools_flagstat' -Cpu 2 -MemoryGb $StageMemory.samtools -OutputPaths @($flagstatHost) }

$countHost = Join-Path $CountDir "$SampleId.featureCounts.txt"
$countSummaryHost = "$countHost.summary"
if (-not (Test-Path -LiteralPath $countHost) -or -not (Test-Path -LiteralPath $countSummaryHost)) {
  $fcArgs = @('featureCounts', '-T', "$Threads", '-a', "/data/reference/$($Config.reference.gtf)",
    '-o', "/data/results/$ResultNamespace/counts/$SampleId.featureCounts.txt", '-s', "$($Config.library.featurecounts_s)")
  if ($Layout -eq 'PAIRED') { $fcArgs += @('-p', '--countReadPairs') }
  $fcArgs += $bamContainer
  Invoke-MeasuredDocker -Stage 'featurecounts_s2' -Image $Images.subread -MemoryGb $StageMemory.featurecounts -Cpu $Threads `
    -LogPath (Join-Path $LogRoot 'featurecounts.log') -Arguments $fcArgs -OutputPaths @($countHost, $countSummaryHost)
} else { Record-Reused -Stage 'featurecounts_s2' -Cpu $Threads -MemoryGb $StageMemory.featurecounts -OutputPaths @($countHost, $countSummaryHost) }

$junctionHost = Join-Path $JunctionDir "$SampleId.junc"
if (-not (Test-Path -LiteralPath $junctionHost) -or -not $junctionPolicyMatches) {
  Invoke-MeasuredDocker -Stage "regtools_$($RegtoolsStrand.ToLowerInvariant())" -Image $Images.regtools -MemoryGb $StageMemory.regtools -Cpu 2 `
    -LogPath (Join-Path $LogRoot 'regtools.log') `
    -Arguments @('regtools', 'junctions', 'extract', '-a', "$($Config.library.regtools_anchor)", '-m', "$($Config.library.regtools_min_intron)", '-M', "$($Config.library.regtools_max_intron)", '-s', $RegtoolsStrand,
      $bamContainer, '-o', "/data/results/$ResultNamespace/junctions/$SampleId.junc") `
    -OutputPaths @($junctionHost)
  $checkArgs = @((Join-Path $PSScriptRoot 'junction_strand_check.py'), '--junctions', $junctionHost,
    '--report', (Join-Path $JunctionDir "$SampleId.strand_check.json"))
  if ($Layout -eq 'SINGLE') { $checkArgs += '--require-both-strands' }
  & python @checkArgs
  if ($LASTEXITCODE -ne 0) { throw "Junction strand check failed for $SampleId; BAM is retained." }
  [ordered]@{ signature = $JunctionSignature; srr = $Srr; layout = $Layout; strand = $RegtoolsStrand;
    junction_sha256 = (Get-FileHash -LiteralPath $junctionHost -Algorithm SHA256).Hash.ToLowerInvariant() } |
    ConvertTo-Json | Set-Content -LiteralPath $JunctionPolicyPath -Encoding utf8
} else { Record-Reused -Stage "regtools_$($RegtoolsStrand.ToLowerInvariant())" -Cpu 2 -MemoryGb $StageMemory.regtools -OutputPaths @($junctionHost) }

$rseqcHost = Join-Path $RseqcDir "$SampleId.infer_experiment.txt"
if (-not (Test-Path -LiteralPath $rseqcHost)) {
  Invoke-MeasuredDocker -Stage 'rseqc_infer_experiment' -Image $Images.rseqc -MemoryGb $StageMemory.rseqc -Cpu 2 `
    -LogPath (Join-Path $LogRoot 'rseqc.log') `
    -Arguments @('bash', '-lc', "infer_experiment.py -r /data/reference/$($Config.library.rseqc_bed) -i '$bamContainer' > '/data/results/$ResultNamespace/rseqc/$SampleId.infer_experiment.txt'") `
    -OutputPaths @($rseqcHost)
} else { Record-Reused -Stage 'rseqc_infer_experiment' -Cpu 2 -MemoryGb $StageMemory.rseqc -OutputPaths @($rseqcHost) }

if ($SampleId -in @($Config.reprocheck.sample_id, $Config.reprocheck_se.sample_id)) {
  $inputDir = Join-Path $ResultRoot 'input_checksums'
  New-Item -ItemType Directory -Path $inputDir -Force | Out-Null
  $inputPath = Join-Path $inputDir "$SampleId.json"
  $hashPaths = @($hostR1)
  if ($Layout -eq 'PAIRED') { $hashPaths += $hostR2 }
  $inputJson = & python (Join-Path $PSScriptRoot 'hash_fastq_content.py') @hashPaths
  if ($LASTEXITCODE -ne 0) { throw 'Primary FASTQ hashing failed; intermediates are retained.' }
  $inputRecords = @((($inputJson | Out-String | ConvertFrom-Json).PSObject.Properties.Value))
  $inputRecord = [ordered]@{ sample_id = $SampleId; srr = $Srr; layout = $Layout; fastq = $inputRecords[0] }
  if ($Layout -eq 'PAIRED') { $inputRecord.fastq_r2 = $inputRecords[1] }
  $inputRecord | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $inputPath -Encoding utf8
}

if (-not $KeepBam) {
  Remove-Item -LiteralPath $bamHost, $baiHost -Force -ErrorAction SilentlyContinue
  foreach ($suffix in @('_STARtmp', '_STARpass1', '_STARgenome')) { Remove-WorkDirectory (Join-Path $StarDir "${SampleId}.$suffix") $StarDir }
}
if (-not $KeepFastq) {
  Remove-WorkDirectory $branchFastqDir $WorkRoot
}

$requiredOutputs = @($starLogHost, $countHost, $countSummaryHost, $junctionHost, $rseqcHost) + $fastqcExpected
$missing = @($requiredOutputs | Where-Object { -not (Test-Path -LiteralPath $_) })
if ($missing.Count -gt 0) { throw "Sample completed with missing outputs: $($missing -join ', ')" }
Write-Status 'done' @{ output_files = $requiredOutputs.Count; bam_retained = [bool]$KeepBam; fastq_source = 'sra'; processing_signature = $ProcessingSignature }
Write-Host "[done] $SampleId; free=$(Get-FreeGb) GB"
