"""Read exported image/PDF/SVG headers only; never renders or modifies a figure."""
import argparse
import csv
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0, 'E:/YFS/workspace/FattyLiver_generation_analysis/.figure_qa_python')
from PIL import Image
import pymupdf as fitz

p = argparse.ArgumentParser()
p.add_argument('run', type=Path)
a = p.parse_args()
with (a.run / 'tables/Figure_export_manifest_EN.csv').open(encoding='utf-8-sig', newline='') as f:
    manifest = list(csv.DictReader(f))
results = []
for item in manifest:
    stem = item['figure']
    w, h = float(item['width_mm']), float(item['height_mm'])
    for suffix in ('png', 'tiff'):
        path = a.run / 'figures' / f'{stem}.{suffix}'
        with Image.open(path) as im:
            dx, dy = im.info.get('dpi', (0, 0))
            ok = (abs(float(dx) - 600) < .02 and abs(float(dy) - 600) < .02
                  and abs(im.width - w / 25.4 * 600) <= 1.1
                  and abs(im.height - h / 25.4 * 600) <= 1.1)
            detail = f'{im.width}x{im.height}px; dpi={dx},{dy}'
            results.append(dict(figure=stem, format=suffix, passed=ok, detail=detail))
    with fitz.open(a.run / 'figures' / f'{stem}.pdf') as doc:
        rect = doc[0].rect
        # R's Cairo PDF device exposes integer-big-point device dimensions
        # (dev.size('in') * 72), truncating the requested physical dimensions.
        # Check that exact device convention; SVG retains fractional points.
        pw, ph = math.floor(w / 25.4 * 72), math.floor(h / 25.4 * 72)
        ok = (len(doc) == 1 and abs(rect.width - pw) < .01
              and abs(rect.height - ph) < .01)
        results.append(dict(figure=stem, format='pdf', passed=ok,
                            detail=f'{rect.width}x{rect.height}pt; Cairo device expected={pw}x{ph}pt; pages={len(doc)}'))
    root = ET.parse(a.run / 'figures' / f'{stem}.svg').getroot()
    sw, sh = root.attrib['width'], root.attrib['height']
    assert sw.endswith('pt') and sh.endswith('pt')
    ok = abs(float(sw[:-2]) - w / 25.4 * 72) < .1 and abs(float(sh[:-2]) - h / 25.4 * 72) < .1
    results.append(dict(figure=stem, format='svg', passed=ok, detail=f'{sw}x{sh}'))
with (a.run / 'audit/export_header_validation.csv').open('w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['figure', 'format', 'passed', 'detail'])
    writer.writeheader()
    writer.writerows(results)
failed = [row for row in results if not row['passed']]
print(f'Export header validation: {len(results)-len(failed)}/{len(results)} PASS')
if failed:
    print(failed)
    raise SystemExit(1)
