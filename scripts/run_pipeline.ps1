param(
  [string]$ConfigPath = '',
  [string]$HostData = '',
  [int]$Threads = 0,
  [double]$MinFreeGb = 0,
  [switch]$RunStrandBenchmark = $true
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'common.ps1')
$ConfigPath = Get-CanonicalConfigPath $ConfigPath
$Config = Get-BranchConfig $ConfigPath
$HostData = Get-DataRoot $Config $HostData
if ($Threads -le 0) { $Threads = [int]$Config.resources.threads }
if ($MinFreeGb -le 0) { $MinFreeGb = [double]$Config.resources.minimum_free_gb }
New-Item -ItemType Directory -Force -Path $HostData | Out-Null

# Each stage verifies its prerequisites before returning.
& (Join-Path $PSScriptRoot 'prepare_reference.ps1') -ConfigPath $ConfigPath -HostData $HostData -Threads $Threads
& (Join-Path $PSScriptRoot 'process_samples.ps1') -ConfigPath $ConfigPath -HostData $HostData -Threads $Threads -MinFreeGb $MinFreeGb
& (Join-Path $PSScriptRoot 'run_repeat.ps1') -Layout PAIRED -ConfigPath $ConfigPath -HostData $HostData -Threads $Threads -RunStrandBenchmark:$RunStrandBenchmark
& (Join-Path $PSScriptRoot 'run_repeat.ps1') -Layout SINGLE -ConfigPath $ConfigPath -HostData $HostData -Threads $Threads
& python (Join-Path $PSScriptRoot 'build_report.py') --branch-root (Get-BranchRoot) --config $ConfigPath --data-root $HostData
if ($LASTEXITCODE -ne 0) { throw 'Report generation failed.' }
