"""Validate the R cohort overview, then use the unique shared results viewer.

No statistical analysis is run and no legacy duplicate viewing tree is created.
"""
from pathlib import Path
import hashlib
import json
import math
import subprocess
import sys
import xml.etree.ElementTree as ET

root = Path(__file__).resolve().parents[1]
run = (root / 'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt').read_text().strip()
current = root / 'outputs/24_Two_generation_reanalysis/current' / run
folder = 'Cohort_overview_two_generation_20261024'
source = current / 'figures' / folder
audit = current / 'audit' / folder
stem = 'Study_population_and_analysis_sets'
sys.path.insert(0, str(root / '.figure_qa_python'))
import pymupdf
from PIL import Image


def digest(path):
    return hashlib.md5(path.read_bytes()).hexdigest()


contract = json.loads((audit / 'figure_contract.json').read_text())
counts = json.loads((audit / 'verified_counts.json').read_text())
width, height = contract['width_mm'], contract['height_mm']
for suffix in ['effective-text-audit.json', 'collision-audit.json']:
    report = json.loads((audit / f'{stem}.{suffix}').read_text())
    assert report['verdict'] == 'PASS', report
with pymupdf.open(source / f'{stem}.pdf') as doc:
    assert len(doc) == 1
    page = doc[0]
    assert page.rect.width == math.floor(width / 25.4 * 72)
    assert page.rect.height == math.floor(height / 25.4 * 72)
    content = page.get_text()
    for number in [counts[k] for k in [
            'clinical_source', 'clinical_ultrasound', 'clinical_cohort',
            'clinical_families', 'microbiome_matched', 'microbiome_QC',
            'microbiome_cohort', 'microbiome_families']
        ] + counts['clinical_generations'] + counts['microbiome_generations']:
        assert f'{number:,}' in content, number
    for forbidden in ['Stool sample not required', 'Minors retained',
                      'Complete cases', 'multiple imputation', 'Model 1',
                      'Model 2', 'Model 3', 'Model 4', 'Model 5', 'Model 6',
                      'P1', 'P2', 'B1', 'B2', '1,275', '319 primary families']:
        assert forbidden not in content, forbidden
    for expected in ['Intergenerational', 'steatosis associations',
                     'Family microbiome', 'similarity',
                     'Steatosis–microbiome', 'associations']:
        assert expected in content, expected
    # One flow schematic: check actual rendered heading alignment, matching
    # cohort geometry, and analysis-box widths rather than claiming panel QA.
    headers = [max(page.search_for(s), key=lambda r: r.y0) for s in [
        'Source data', 'Study populations', 'Analyses']]
    assert max(r.y0 for r in headers) - min(r.y0 for r in headers) < .01
    drawings = page.get_drawings()
    cohort_boxes = [d['rect'] for d in drawings
                    if abs(d['rect'].width - 52 / 25.4 * 72) < .1
                    and abs(d['rect'].height - 43 / 25.4 * 72) < .1]
    analysis_boxes = [d['rect'] for d in drawings
                      if abs(d['rect'].width - 53 / 25.4 * 72) < .1]
    assert len(cohort_boxes) == 2 and len(analysis_boxes) == 3
    assert abs(cohort_boxes[0].x0 - cohort_boxes[1].x0) < .01
    assert max(r.x0 for r in analysis_boxes) - min(r.x0 for r in analysis_boxes) < .01
    geometry = {
        'heading_baseline_spread_pt': max(r.y0 for r in headers) - min(r.y0 for r in headers),
        'cohort_boxes_pt': [list(r) for r in cohort_boxes],
        'analysis_boxes_pt': [list(r) for r in analysis_boxes],
        'classification': 'One unified flow schematic; quantitative panel alignment not applicable',
        'verdict': 'PASS'}
(audit / 'rendered_geometry.json').write_text(json.dumps(geometry, indent=2))
for ext in ['png', 'tiff']:
    with Image.open(source / f'{stem}.{ext}') as im:
        assert abs(im.width - width / 25.4 * 600) < 2
        assert abs(im.height - height / 25.4 * 600) < 2
        assert min(im.info['dpi']) > 599
ET.parse(source / f'{stem}.svg')

title = 'Study populations and analyses'
(source / 'index.html').write_text(f'''<!doctype html><meta charset="utf-8"><title>{title}</title>
<style>body{{font:16px/1.5 Arial,sans-serif;max-width:1150px;margin:25px auto;color:#171717}}img{{width:100%}}a{{color:#245882}}</style>
<h1>{title}</h1>
<p><a href="{stem}.pdf">PDF</a> | <a href="{stem}.png">PNG</a> | <a href="{stem}.tiff">TIFF</a> | <a href="{stem}.svg">SVG</a> | <a href="Source_counts_EN.csv">Source counts</a> | <a href="Figure_legend_EN.txt">English legend</a></p>
<img src="{stem}.png" alt="Source data, clinical and microbiome cohorts, and three analysis directions">
''', encoding='utf-8')

# Use the study-ordered shared viewer, never the retired duplicate directory.
subprocess.run([str(root / 'tools/R-4.5.2/bin/x64/Rscript.exe'),
                str(root / 'scripts/refresh_results_gallery.R')], cwd=root, check=True)
dest = (root / 'Results_表格与图片/01_当前分析/01_基础队列与临床及群落分析/figures' / folder)
for f in source.iterdir():
    if f.is_file():
        assert digest(dest / f.name) == digest(f), f.name
report = {'status': 'PASS', 'formats': ['PDF', 'PNG', 'TIFF', 'SVG'],
          'clinical_N': counts['clinical_cohort'], 'microbiome_N': counts['microbiome_cohort'],
          'output': str(dest), 'no_statistics_rerun': True,
          'change': 'Both marked cohort notes removed; right column reduced to three established analysis directions.',
          'visual_review': 'Final R-rendered image inspected for readable labels, correct counts, aligned boxes, and clear connectors.'}
(audit / 'delivery_validation.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2, ensure_ascii=False))
