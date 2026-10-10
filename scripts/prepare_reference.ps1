param(
  [string]$ConfigPath = '',
  [string]$HostData = '',
  [int]$Threads = 0,
  [int]$StarMemoryGb = 0
)

$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'common.ps1')
$ConfigPath = Get-CanonicalConfigPath $ConfigPath
$Config = Get-BranchConfig $ConfigPath
$HostData = Get-DataRoot $Config $HostData
if ($Threads -le 0) { $Threads = [int]$Config.resources.threads }
if ($StarMemoryGb -le 0) { $StarMemoryGb = [int]$Config.resources.star_memory_gb }
$ReferenceDir = Join-Path $HostData 'reference'
$LogDir = Join-Path (Get-DataPath $HostData $Config.paths.data_logs) 'reference'
$StatusPath = Join-Path $LogDir 'status.json'
$GtfGz = Join-Path $ReferenceDir $Config.reference.gtf_gz
$Gtf = Join-Path $ReferenceDir $Config.reference.gtf
$Fasta = Join-Path $ReferenceDir $Config.reference.fasta
$FinalIndex = Join-Path $ReferenceDir $Config.reference.star_index
$BuildIndex = Join-Path $ReferenceDir ($Config.reference.star_index + '_building')
$Url = $Config.reference.gtf_url
$Md5Url = $Config.reference.md5sums_url
$Md5Sums = Join-Path $ReferenceDir 'gencode.v50.MD5SUMS'
$ExpectedGtfSha256 = $Config.reference.gtf_sha256
$PigzImage = $Config.containers.pigz
$StarImage = $Config.containers.star
$referenceBoundary = [IO.Path]::GetFullPath($ReferenceDir).TrimEnd('\') + '\'
foreach ($path in @($GtfGz, $Gtf, $Fasta, $BuildIndex, $FinalIndex)) {
  if (-not [IO.Path]::GetFullPath($path).StartsWith($referenceBoundary, [StringComparison]::OrdinalIgnoreCase)) { throw 'Reference target is outside the reference directory.' }
}

New-Item -ItemType Directory -Force -Path $ReferenceDir, $LogDir | Out-Null

function Write-Status([string]$Phase, [hashtable]$Extra = @{}) {
  $payload = @{
    phase = $Phase
    time = (Get-Date).ToString('o')
    free_gb = Get-FreeSpaceGb $HostData
  }
  foreach ($key in $Extra.Keys) { $payload[$key] = $Extra[$key] }
  $payload | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $StatusPath -Encoding utf8
}

if (-not (Test-Path -LiteralPath $Fasta)) {
  Write-Status 'download_fasta'
  $fastaGz = Join-Path $ReferenceDir $Config.reference.fasta_gz
  & curl.exe --location --fail --retry 5 --continue-at - --output $fastaGz $Config.reference.fasta_url
  if ($LASTEXITCODE -ne 0) { throw 'FASTA download failed.' }
  $fastaBuilding = "$Fasta.building"
  $cmd = "pigz -dc /data/reference/$($Config.reference.fasta_gz) > /data/reference/$($Config.reference.fasta).building"
  & docker run --rm --memory="$($Config.resources.stage_memory_gb.compression)g" --cpus=2 -v "${HostData}:/data" $PigzImage bash -lc $cmd *> (Join-Path $LogDir 'decompress_fasta.log')
  if ($LASTEXITCODE -ne 0) { throw 'FASTA decompression failed.' }
  $hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $fastaBuilding).Hash.ToLowerInvariant()
  if ($hash -ne $Config.reference.fasta_sha256) { throw 'FASTA SHA-256 mismatch.' }
  Move-Item -LiteralPath $fastaBuilding -Destination $Fasta
}
$FastaHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $Fasta).Hash.ToLowerInvariant()
if ($FastaHash -ne $Config.reference.fasta_sha256) { throw 'Frozen FASTA SHA-256 mismatch.' }

Write-Status 'download_gtf'
if (-not (Test-Path -LiteralPath $GtfGz) -or (Get-FileHash -Algorithm MD5 -LiteralPath $GtfGz).Hash.ToLowerInvariant() -ne $Config.reference.gtf_gz_md5) {
  & curl.exe --location --fail --retry 5 --continue-at - --output $GtfGz $Url
  if ($LASTEXITCODE -ne 0) { throw "GTF download failed with exit $LASTEXITCODE" }
}

$downloadedGtfMd5 = (Get-FileHash -Algorithm MD5 -LiteralPath $GtfGz).Hash.ToLowerInvariant()
if ($downloadedGtfMd5 -ne $Config.reference.gtf_gz_md5) { throw 'Downloaded GTF checksum mismatch.' }
Write-Status 'decompress_gtf'
if (-not (Test-Path -LiteralPath $Gtf) -or (Get-Item -LiteralPath $Gtf).Length -lt 1GB) {
  $decompressLog = Join-Path $LogDir 'decompress.log'
  $gtfBuilding = "$Gtf.building"
  Remove-Item -LiteralPath $gtfBuilding -Force -ErrorAction SilentlyContinue
  $cmd = "pigz -dc /data/reference/$($Config.reference.gtf_gz) > /data/reference/$($Config.reference.gtf).building"
  & docker run --rm --memory="$($Config.resources.stage_memory_gb.compression)g" --cpus=2 -v "${HostData}:/data" $PigzImage bash -lc $cmd *> $decompressLog
  if ($LASTEXITCODE -ne 0) { throw "GTF decompression failed; see $decompressLog" }
  Move-Item -LiteralPath $gtfBuilding -Destination $Gtf
}

$GtfHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $Gtf).Hash.ToLowerInvariant()
if ($GtfHash -ne $ExpectedGtfSha256) { throw "Uncompressed GTF SHA-256 mismatch: observed $GtfHash, expected $ExpectedGtfSha256" }
$GtfGzHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $GtfGz).Hash.ToLowerInvariant()
$GtfGzMd5 = (Get-FileHash -Algorithm MD5 -LiteralPath $GtfGz).Hash.ToLowerInvariant()
if (-not (Test-Path -LiteralPath $Md5Sums)) {
  & curl.exe --location --fail --retry 5 --output $Md5Sums $Md5Url
  if ($LASTEXITCODE -ne 0) { throw "GENCODE MD5SUMS download failed with exit $LASTEXITCODE" }
}
$publishedLine = Get-Content -LiteralPath $Md5Sums | Where-Object { $_ -match '\s+gencode\.v50\.primary_assembly\.annotation\.gtf\.gz$' } | Select-Object -First 1
if (-not $publishedLine) { throw 'Final GTF is absent from the published GENCODE MD5SUMS file' }
$PublishedMd5 = ($publishedLine -split '\s+')[0].ToLowerInvariant()
if ($GtfGzMd5 -ne $PublishedMd5) { throw "GTF MD5 mismatch: observed $GtfGzMd5, published $PublishedMd5" }
if ($GtfGzMd5 -ne $Config.reference.gtf_gz_md5) { throw "GTF MD5 differs from frozen config: observed $GtfGzMd5, configured $($Config.reference.gtf_gz_md5)" }
@(
  "$GtfHash  $($Config.reference.gtf)",
  "$GtfGzHash  $($Config.reference.gtf_gz)"
) | Set-Content -LiteralPath (Join-Path $ReferenceDir 'BRANCH_A_GTF_SHA256SUMS') -Encoding ascii

$required = @('Genome', 'SA', 'SAindex', 'genomeParameters.txt')
$ready = $true
foreach ($name in $required) {
  if (-not (Test-Path -LiteralPath (Join-Path $FinalIndex $name))) { $ready = $false }
}
if ($ready) {
  $parameters = Get-Content -LiteralPath (Join-Path $FinalIndex 'genomeParameters.txt') -Raw
  foreach ($pair in @(@('sjdbOverhang', $Config.reference.sjdbOverhang), @('genomeSAsparseD', $Config.reference.genomeSAsparseD))) {
    if ($parameters -notmatch ('(?m)^' + $pair[0] + '\s+' + $pair[1] + '\s*$')) { throw 'Existing STAR index parameters differ from config.' }
  }
  $manifestPath = Join-Path $ReferenceDir 'BRANCH_A_REFERENCE_READY.json'
  if (Test-Path -LiteralPath $manifestPath) {
    $existingManifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
    $existingManifest | Add-Member -NotePropertyName gtf_gz_md5 -NotePropertyValue $GtfGzMd5 -Force
    $existingManifest | Add-Member -NotePropertyName published_gtf_gz_md5 -NotePropertyValue $PublishedMd5 -Force
    $existingManifest | Add-Member -NotePropertyName md5_verified -NotePropertyValue $true -Force
    $existingManifest | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $manifestPath -Encoding utf8
  } else { throw 'Existing index lacks its reference manifest; choose a new index directory.' }
  if ($existingManifest.fasta_sha256 -ne $FastaHash -or $existingManifest.gtf_sha256 -ne $GtfHash) { throw 'Existing index reference hashes differ from config.' }
  Write-Status 'done' @{ gtf_sha256 = $GtfHash; index = $FinalIndex; reused = $true }
  Write-Host "Final methods reference already ready: $FinalIndex"
  exit 0
}

if (Test-Path -LiteralPath $BuildIndex) {
  $active = & docker ps --filter 'name=^/altersplcinf-a-reference$' --format '{{.Names}}'
  if ($LASTEXITCODE -ne 0) { throw 'Docker is unavailable; the incomplete index is retained.' }
  if ($active) { throw 'Reference preparation is already running.' }
  Remove-WorkDirectory $BuildIndex $ReferenceDir
}
New-Item -ItemType Directory -Path $BuildIndex | Out-Null

Write-Status 'build_star_index' @{ gtf_sha256 = $GtfHash }
$buildLog = Join-Path $LogDir 'star_genomeGenerate.log'
$started = Get-Date
& docker run --rm --name altersplcinf-a-reference --memory="${StarMemoryGb}g" --cpus=$Threads `
  -v "${HostData}:/data" -w /data/reference $StarImage STAR `
  --runMode genomeGenerate `
  --runThreadN $Threads `
  --genomeDir "/data/reference/$($Config.reference.star_index)_building" `
  --genomeFastaFiles "/data/reference/$($Config.reference.fasta)" `
  --sjdbGTFfile "/data/reference/$($Config.reference.gtf)" `
  --sjdbOverhang $Config.reference.sjdbOverhang `
  --genomeSAsparseD $Config.reference.genomeSAsparseD `
  --limitGenomeGenerateRAM $Config.reference.limit_genome_generate_ram_bytes *> $buildLog
$exitCode = $LASTEXITCODE
$elapsed = [math]::Round(((Get-Date) - $started).TotalSeconds, 1)
if ($exitCode -ne 0) {
  Write-Status 'failed' @{ step = 'build_star_index'; exit_code = $exitCode; elapsed_seconds = $elapsed }
  throw "STAR genomeGenerate failed; see $buildLog"
}

foreach ($name in $required) {
  if (-not (Test-Path -LiteralPath (Join-Path $BuildIndex $name))) {
    throw "STAR build completed without required file: $name"
  }
}
Move-Item -LiteralPath $BuildIndex -Destination $FinalIndex

$indexBytes = (Get-ChildItem -LiteralPath $FinalIndex -File | Measure-Object Length -Sum).Sum
$manifest = @{
  ready = $true
  completed = (Get-Date).ToString('o')
  fasta = $Config.reference.fasta
  fasta_sha256 = (Get-FileHash -Algorithm SHA256 -LiteralPath $Fasta).Hash.ToLowerInvariant()
  gtf = $Config.reference.gtf
  gtf_sha256 = $GtfHash
  gtf_gz_sha256 = $GtfGzHash
  gtf_gz_md5 = $GtfGzMd5
  published_gtf_gz_md5 = $PublishedMd5
  md5_verified = $true
  star_index = $Config.reference.star_index
  sjdbOverhang = $Config.reference.sjdbOverhang
  genomeSAsparseD = $Config.reference.genomeSAsparseD
  limit_genome_generate_ram_bytes = $Config.reference.limit_genome_generate_ram_bytes
  threads = $Threads
  memory_limit_gb = $StarMemoryGb
  elapsed_seconds = $elapsed
  index_size_gb = [math]::Round($indexBytes / 1GB, 3)
}
$manifest | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $ReferenceDir 'BRANCH_A_REFERENCE_READY.json') -Encoding utf8
Write-Status 'done' @{ gtf_sha256 = $GtfHash; index = $FinalIndex; elapsed_seconds = $elapsed; reused = $false }
Write-Host "Methods reference ready: $FinalIndex"
