# One-time checkpoint handoff for this authorized long SGB computation.
# Explicit PID belongs to the Rscript process started in this task.
param([int]$SgbProcessId=4008,[string]$RunId='sgb_associations_v1_20261022')
$ErrorActionPreference='Stop'
$sgbRoot=Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $sgbRoot
if($RunId -notmatch '^[A-Za-z0-9][A-Za-z0-9_.-]*$'){throw 'Unsafe run ID'}
$ledger=Join-Path $sgbRoot "outputs/17_SGB_association/_running/$RunId/audit/overall_progress.csv"
$sgbR=Join-Path $sgbRoot 'tools/R-4.5.2/bin/Rscript.exe'
$expectedBin=Join-Path $sgbRoot 'tools/R-4.5.2/bin'
while(-not(Test-Path -LiteralPath $ledger)){
  if(-not(Get-Process -Id $SgbProcessId -ErrorAction SilentlyContinue)){throw 'Original process exited before checkpoint'}
  Start-Sleep -Seconds 2
}
$done=Import-Csv -LiteralPath $ledger
if(-not($done|Where-Object {$_.status -eq 'COMPLETE'})){throw 'No completed batch in ledger'}
$running=Get-CimInstance Win32_Process
$main=$running|Where-Object {$_.ProcessId -eq $SgbProcessId}
if(-not $main -or -not $main.ExecutablePath.StartsWith($expectedBin,[StringComparison]::OrdinalIgnoreCase) -or $main.CommandLine -notmatch '17_sgb_association\.R'){
  throw 'Refusing to stop a process outside this task'
}
$ids=@([uint32]$SgbProcessId)
do{
  $new=@($running|Where-Object {$_.ParentProcessId -in $ids -and $_.ProcessId -notin $ids}|Select-Object -ExpandProperty ProcessId)
  $ids+= $new
}while($new.Count -gt 0)
# Freeze the parent first so it cannot dispatch more work, then stop its known
# R worker descendants. Completed signed batch files are not removed.
Stop-Process -Id $SgbProcessId -Force
foreach($idItem in ($ids|Where-Object {$_ -ne $SgbProcessId})){
  $child=$running|Where-Object {$_.ProcessId -eq $idItem}
  if($child.ExecutablePath -and $child.ExecutablePath.StartsWith($expectedBin,[StringComparison]::OrdinalIgnoreCase)){
    Stop-Process -Id $idItem -Force -ErrorAction SilentlyContinue
  }
}
Add-Content -LiteralPath 'tools/downloads/sgb_checkpoint_handoff.log' -Value "$(Get-Date -Format o) resumed after a complete signed batch; entry/code updated, same statistical models."
$env:FATTYLIVER_SGB_RESUME_RUN=$RunId
& $sgbR 'scripts/17_sgb_association.R' '--run' *> 'tools/downloads/sgb_resumed_run.log'
if($LASTEXITCODE -ne 0){throw "Resumed SGB run failed: $LASTEXITCODE; see sgb_resumed_run.log"}
