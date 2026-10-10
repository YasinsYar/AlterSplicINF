function Get-BranchConfig {
  param([string]$ConfigPath)
  $loader = Join-Path $PSScriptRoot 'load_config.py'
  $json = & python $loader --config $ConfigPath
  if ($LASTEXITCODE -ne 0 -or -not $json) { throw "Unable to load branch config: $ConfigPath" }
  $config = $json | Out-String | ConvertFrom-Json
  if ($config.library.featurecounts_s -ne 2 -or $config.library.regtools_strand -ne 'RF' -or
      $config.library.regtools_strand_se -ne 'XS' -or $config.library.star_twopass_mode -ne 'Basic' -or
      $config.library.star_out_sam_strand_field -ne 'intronMotif') {
    throw 'The published protocol requires featureCounts s2, PE RF, SE XS and STAR Basic/intronMotif.'
  }
  return $config
}

function Get-CanonicalConfigPath {
  param([string]$RequestedPath)
  if ($RequestedPath) { return (Resolve-Path -LiteralPath $RequestedPath).Path }
  return (Resolve-Path -LiteralPath (Join-Path (Split-Path -Parent $PSScriptRoot) 'config.yaml')).Path
}

function Get-BranchRoot {
  return [IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
}

function Get-DataRoot {
  param($Config, [string]$RequestedPath = '')
  if ($RequestedPath) { return [IO.Path]::GetFullPath($RequestedPath) }
  $configured = [string]$Config.paths.data_root
  if ([IO.Path]::IsPathRooted($configured)) { return [IO.Path]::GetFullPath($configured) }
  return [IO.Path]::GetFullPath((Join-Path (Get-BranchRoot) $configured))
}

function Get-DataPath {
  param([string]$DataRoot, [string]$ConfiguredPath)
  if ([IO.Path]::IsPathRooted($ConfiguredPath)) { return [IO.Path]::GetFullPath($ConfiguredPath) }
  return [IO.Path]::GetFullPath((Join-Path $DataRoot $ConfiguredPath))
}

function Get-FreeSpaceGb {
  param([string]$Path)
  $root = [IO.Path]::GetPathRoot([IO.Path]::GetFullPath($Path))
  $driveName = $root.TrimEnd([IO.Path]::DirectorySeparatorChar, [IO.Path]::AltDirectorySeparatorChar).TrimEnd(':')
  return [math]::Round((Get-PSDrive -Name $driveName).Free / 1GB, 3)
}

function Get-RegtoolsStrand {
  param($Config, [ValidateSet('PAIRED', 'SINGLE')][string]$Layout)
  if ($Layout -eq 'SINGLE') {
    if ($Config.library.regtools_strand_se -ne 'XS' -or $Config.library.star_out_sam_strand_field -ne 'intronMotif') {
      throw 'Single-end junction extraction requires XS mode and STAR intronMotif strand tags.'
    }
    return [string]$Config.library.regtools_strand_se
  }
  return [string]$Config.library.regtools_strand
}

function Get-JunctionPolicySignature {
  param($Config, [ValidateSet('PAIRED', 'SINGLE')][string]$Layout)
  $policy = [ordered]@{
    layout = $Layout
    strand = Get-RegtoolsStrand $Config $Layout
    anchor = $Config.library.regtools_anchor
    min_intron = $Config.library.regtools_min_intron
    max_intron = $Config.library.regtools_max_intron
    regtools_image = $Config.containers.regtools
    star_image = $Config.containers.star
    star_mode = $Config.library.star_twopass_mode
    star_strand_field = $Config.library.star_out_sam_strand_field
    fasta_sha256 = $Config.reference.fasta_sha256
    gtf_sha256 = $Config.reference.gtf_sha256
  }
  $bytes = [Text.Encoding]::UTF8.GetBytes(($policy | ConvertTo-Json -Compress))
  $digest = [Security.Cryptography.SHA256]::Create()
  try { return ([BitConverter]::ToString($digest.ComputeHash($bytes))).Replace('-', '').ToLowerInvariant() }
  finally { $digest.Dispose() }
}

function Test-CompleteBam {
  param([string]$Path)
  if (-not (Test-Path -LiteralPath $Path -PathType Leaf) -or (Get-Item -LiteralPath $Path).Length -le 28) { return $false }
  # A completed BGZF BAM ends with the standard 28-byte EOF block.
  $expected = [byte[]]@(31,139,8,4,0,0,0,0,0,255,6,0,66,67,2,0,27,0,3,0,0,0,0,0,0,0,0,0)
  $stream = [IO.File]::OpenRead($Path)
  try {
    $tail = [byte[]]::new(28)
    $stream.Seek(-28, [IO.SeekOrigin]::End) | Out-Null
    return $stream.Read($tail, 0, 28) -eq 28 -and [Convert]::ToBase64String($tail) -eq [Convert]::ToBase64String($expected)
  } finally { $stream.Dispose() }
}

function Remove-WorkDirectory {
  param([string]$Path, [string]$Root)
  $target = [IO.Path]::GetFullPath($Path).TrimEnd('\')
  $boundary = [IO.Path]::GetFullPath($Root).TrimEnd('\') + '\'
  if (-not $target.StartsWith($boundary, [StringComparison]::OrdinalIgnoreCase)) {
    throw "Cleanup target is outside the work directory: $target"
  }
  if (Test-Path -LiteralPath $target) { Remove-Item -LiteralPath $target -Recurse -Force }
}

function Get-ProcessingSignature {
  param($Config, [string]$Srr, [string]$Layout)
  $policy = [ordered]@{ srr = $Srr; layout = $Layout; reference = $Config.reference; library = $Config.library; containers = $Config.containers }
  $bytes = [Text.Encoding]::UTF8.GetBytes(($policy | ConvertTo-Json -Depth 8 -Compress))
  $digest = [Security.Cryptography.SHA256]::Create()
  try { return ([BitConverter]::ToString($digest.ComputeHash($bytes))).Replace('-', '').ToLowerInvariant() }
  finally { $digest.Dispose() }
}
