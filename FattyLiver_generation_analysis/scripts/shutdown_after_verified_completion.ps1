param([Parameter(Mandatory=$true)][string]$CompletionPath)
$ErrorActionPreference = 'Stop'
$projectPath = 'E:\YFS\workspace\FattyLiver_generation_analysis'
$projectResolved = (Resolve-Path -LiteralPath $projectPath).Path
$completionResolved = (Resolve-Path -LiteralPath $CompletionPath).Path
if (-not $completionResolved.StartsWith($projectResolved + [IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)) {
    throw 'Completion manifest must be inside the fatty-liver project.'
}
$record = Get-Content -LiteralPath $completionResolved -Raw -Encoding UTF8 | ConvertFrom-Json
if ($record.status -ne 'COMPLETED' -or $record.shutdown_authorized -ne $true -or $record.validation_passed -ne $true -or $record.figure_qa_passed -ne $true) {
    throw 'Work, statistical validation or figure QA is incomplete; shutdown was not requested.'
}
$logPath = Join-Path $projectResolved 'shutdown_approved_design.log'
Add-Content -LiteralPath $logPath -Encoding UTF8 -Value ('Verified completion; non-forced shutdown helper started at ' + (Get-Date -Format o))
# User explicitly requested shutdown after completion. No /f is used; other applications can block shutdown.
Start-Sleep -Seconds 60
if (Test-Path -LiteralPath (Join-Path $projectResolved 'CANCEL_SHUTDOWN')) {
    Add-Content -LiteralPath $logPath -Encoding UTF8 -Value 'Cancelled by project CANCEL_SHUTDOWN marker.'
    exit 0
}
foreach ($runPath in $record.completed_runs) {
    $runResolved = (Resolve-Path -LiteralPath ([string]$runPath)).Path
    if (-not $runResolved.StartsWith($projectResolved + [IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)) { throw 'Run path outside project.' }
    $statusPath = Join-Path $runResolved 'run_status.txt'
    if ((Get-Content -LiteralPath $statusPath -Encoding UTF8 -TotalCount 1) -ne 'COMPLETED') { throw 'A required run is incomplete.' }
}
Add-Content -LiteralPath $logPath -Encoding UTF8 -Value ('Invoking shutdown /s /t 0 without forced application closure at ' + (Get-Date -Format o))
& "$env:SystemRoot\System32\shutdown.exe" /s /t 0
if ($LASTEXITCODE -ne 0) { throw ('Shutdown request failed: ' + $LASTEXITCODE) }
