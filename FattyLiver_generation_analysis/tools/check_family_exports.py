"""Read-only final-file QA. No statistics, graphics device or image rendering."""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '.figure_qa_python'))
from PIL import Image
import pymupdf

root = Path(__file__).resolve().parents[1]
reports = []
for module in ('45_Family_de_novo_discovery', '48_Family_functional_explanation'):
    base = root / 'outputs' / module
    run = base / 'current' / (base / 'LATEST_RUN.txt').read_text(encoding='utf-8-sig').strip()
    assert (run / 'run_status.txt').read_text(encoding='utf-8-sig').splitlines()[0] == 'COMPLETED'
    for pdf in sorted((run / 'figures').glob('*.pdf')):
        with pymupdf.open(pdf) as doc:
            assert len(doc) == 1 and len(doc[0].get_text().strip()) > 10
            pdf_points = (doc[0].rect.width, doc[0].rect.height)
        assert pdf.with_suffix('.svg').exists()
        svg = pdf.with_suffix('.svg').read_text(encoding='utf-8')
        assert '<text' in svg, 'Editable SVG text missing'
        element = ET.fromstring(svg)
        assert element.attrib['width'].endswith('pt') and element.attrib['height'].endswith('pt')
        svg_points = tuple(float(element.attrib[n][:-2]) for n in ('width', 'height'))
        # R/Cairo writes integer-point PDF MediaBoxes. Check this known <1 pt
        # quantization explicitly; validate raster dimensions against the SVG's
        # fractional-point physical contract, rather than the rounded PDF box.
        assert max(abs(pdf_points[i] - svg_points[i]) for i in range(2)) < 1
        width_mm, height_mm = (v * 25.4 / 72 for v in svg_points)
        for ext in ('.png', '.tiff'):
            with Image.open(pdf.with_suffix(ext)) as im:
                dpi = im.info.get('dpi')
                assert dpi and all(abs(float(v) - 600) < 1 for v in dpi), (pdf.name, ext, dpi)
                dpi = tuple(float(v) for v in dpi)
                expected = (width_mm / 25.4 * 600, height_mm / 25.4 * 600)
                assert max(abs(im.size[i] - expected[i]) for i in range(2)) < 2
                reports.append(dict(module=module, figure=pdf.stem, extension=ext,
                                    width_mm=width_mm, height_mm=height_mm,
                                    pdf_media_box_points=pdf_points, svg_points=svg_points,
                                    cairo_pdf_integer_point_quantization_verified=True,
                                    pixels=im.size, dpi=dpi, mode=im.mode,
                                    editable_svg_text=True, original_files_read_only=True))
output = root / 'outputs/family_four_stage_20261104/export_metadata_verification.json'
output.write_text(json.dumps(dict(status='PASS', records=reports), indent=2), encoding='utf-8')
print(f'PASS: {len(reports)} raster files at 600 dpi; PDF dimensions and editable SVG text verified.')
