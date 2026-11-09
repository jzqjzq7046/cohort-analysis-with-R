# Remove only explicitly superseded 25/26/27 implementations after both revised
# authoritative modules have completed. No source-data or module-24 path is used.
$ErrorActionPreference = 'Stop'
$taskProject = [IO.Path]::GetFullPath('E:/YFS/workspace/FattyLiver_generation_analysis').TrimEnd('\','/')
$taskExpected = @{
  '25_Taxa_robustness' = 'taxa_effect_comparisons_v2_20261026'
  '26_Taxon_family_sharing' = 'taxon_family_regression_v2_20261026'
}
foreach ($taskModule in $taskExpected.Keys) {
  $taskId = (Get-Content -LiteralPath "$taskProject/outputs/$taskModule/LATEST_RUN.txt" -TotalCount 1).Trim()
  if ($taskId -ne $taskExpected[$taskModule]) { throw "Revised run is not adopted: $taskModule" }
  $taskStatus = Get-Content -LiteralPath "$taskProject/outputs/$taskModule/current/$taskId/run_status.txt" -TotalCount 1
  if ($taskStatus -ne 'COMPLETED') { throw "Revised run is not complete: $taskModule" }
}
$taskRelativeTargets = @(
  'R/taxa_followup_20261026.R',
  'R/taxa_visualization_scheme_A.R',
  'scripts/25_recover_figures.R',
  'scripts/complete_revised_taxa_figures.R',
  'scripts/27_visualize_taxa_scheme_A.R',
  'scripts/27_redraw_figure25_without_c.R',
  'scripts/27_README_zh_en.md',
  'scripts/taxon_sharing_kernel.py',
  'scripts/__pycache__/taxon_sharing_kernel.cpython-312.pyc',
  'tests/taxa_followup_contract_test.R',
  'tests/test_taxon_sharing_kernel.py',
  'outputs/25_Taxa_robustness/history',
  'outputs/25_Taxa_robustness/_checkpoints/private/9a1e50b63073bdbdf85ea278',
  'outputs/26_Taxon_family_sharing/history',
  'outputs/26_Taxon_family_sharing/_checkpoints/private/b9f79f1f976897e9ac5ed3ee',
  'outputs/27_Taxa_scheme_A_visualization',
  'Results_表格与图片/Followup_taxa_20261026'
)
foreach ($taskCategory in @('Figures_图片','Tables_表格','Reports_报告')) {
  foreach ($taskModule in @('25_Taxa_robustness','26_Taxon_family_sharing','27_Taxa_scheme_A_visualization')) {
    $taskRelativeTargets += "Results_表格与图片/$taskCategory/$taskModule"
  }
}
$taskAudit = @()
foreach ($taskRelative in $taskRelativeTargets) {
  $taskTarget = [IO.Path]::GetFullPath((Join-Path $taskProject $taskRelative))
  if (-not $taskTarget.StartsWith($taskProject + '\', [StringComparison]::OrdinalIgnoreCase)) { throw "Unsafe target: $taskTarget" }
  if (Test-Path -LiteralPath $taskTarget) {
    $taskResolved = (Resolve-Path -LiteralPath $taskTarget).Path
    if (-not $taskResolved.StartsWith($taskProject + '\', [StringComparison]::OrdinalIgnoreCase)) { throw "Resolved path outside project: $taskResolved" }
    Remove-Item -LiteralPath $taskResolved -Recurse -Force
    $taskAudit += [PSCustomObject]@{ relative_path=$taskRelative; removed=$true }
  }
}
# Drop manifest rows for exactly the retired managed copies, so an open Explorer
# folder can be refreshed in place without preserving obsolete generated files.
$taskManifest = "$taskProject/Results_表格与图片/file_manifest.csv"
$taskRows = Import-Csv -LiteralPath $taskManifest
$taskRows | Where-Object { $_.module -notin @('25_Taxa_robustness','26_Taxon_family_sharing','27_Taxa_scheme_A_visualization') } | Export-Csv -LiteralPath $taskManifest -NoTypeInformation -Encoding utf8BOM
$taskAudit | Export-Csv -LiteralPath "$taskProject/outputs/26_Taxon_family_sharing/current/taxon_family_regression_v2_20261026/audit/retired_superseded_paths.csv" -NoTypeInformation -Encoding utf8BOM
Write-Output "Retired $($taskAudit.Count) superseded paths. Revised current runs and source data retained."
