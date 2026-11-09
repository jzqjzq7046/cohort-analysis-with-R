"""Publish the already computed approved-design summary; no refits."""
from pathlib import Path
import shutil
ROOT = Path(__file__).resolve().parents[1]
def main():
    runs = list((ROOT / 'outputs/10_Results_visualization/current').iterdir())
    if len(runs) != 1 or 'approved_results_summary_v1' not in runs[0].name:
        raise RuntimeError('Run scripts/10_visualize_results.R for the approved summary first')
    if (runs[0] / 'run_status.txt').read_text(encoding='utf-8-sig').splitlines()[0] != 'COMPLETED':
        raise RuntimeError('Summary not complete')
    source = runs[0] / 'reports/results_summary_zh.md'
    shutil.copy2(source, ROOT / 'outputs/RESULTS_zh.md')
    print('Approved numerical summary copied from the adopted result; no models refitted')
if __name__ == '__main__':
    main()
