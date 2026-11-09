"""Read-only validation of the completed analysis and its viewing copies."""
from pathlib import Path
import csv, hashlib, json, math, sys, zipfile, xml.etree.ElementTree as ET

root = Path(__file__).resolve().parents[1]
run = (root/'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt').read_text().strip()
source = root/'outputs/24_Two_generation_reanalysis/current'/run
viewer = root/'Results_表格与图片/Reanalysis_two_generation_20261024'/run
assert (source/'run_status.txt').read_text().splitlines()[0] == 'COMPLETED'

def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

figures = rows(source/'tables/Figure_export_manifest_EN.csv')
assert len(figures) == 12 and len({r['figure'] for r in figures}) == 12
sys.path.insert(0, str(root/'.figure_qa_python'))
from PIL import Image
import pymupdf
geometry = []
for r in figures:
    stem = r['figure']
    for ext in ['pdf','png','svg','tiff']:
        assert (source/'figures'/f'{stem}.{ext}').stat().st_size > 100
    with pymupdf.open(source/'figures'/f'{stem}.pdf') as doc:
        assert len(doc) == 1
        # R's Cairo device quantizes MediaBox dimensions down to whole points.
        # Check that exact device rule, rather than requiring fractional points
        # that this backend cannot write. Rendered text/layout QA uses real PDFs.
        expected_width=math.floor(float(r['width_mm'])/25.4*72)
        expected_height=math.floor(float(r['height_mm'])/25.4*72)
        assert abs(doc[0].rect.width-expected_width) < .001
        assert abs(doc[0].rect.height-expected_height) < .001
        geometry.append({'figure':stem,'format':'pdf','MediaBox_pt':list(doc[0].rect),'device_rule':'Cairo whole-point MediaBox; requested mm dimensions recorded in manifest'})
    for ext in ['png','tiff']:
        with Image.open(source/'figures'/f'{stem}.{ext}') as im:
            dpi = im.info.get('dpi')
            assert dpi and min(map(float,dpi)) > 599
            assert abs(im.width - float(r['width_mm'])/25.4*600) < 2
            assert abs(im.height - float(r['height_mm'])/25.4*600) < 2
            geometry.append({'figure':stem,'format':ext,'pixels':list(im.size),'dpi':list(map(float,dpi))})
    ET.parse(source/'figures'/f'{stem}.svg')

for r in rows(viewer/'file_manifest.csv'):
    assert not r['file'].startswith(('private/','audit/'))
    a=source/r['file'];b=viewer/r['file']
    md5=lambda f:hashlib.md5(f.read_bytes()).hexdigest()
    assert md5(a)==r['source_md5']==md5(b)

book=source/'tables/Hepatic_steatosis_two_generation_tables_EN.xlsx'
with zipfile.ZipFile(book) as z:
    assert z.testzip() is None
    doc=ET.fromstring(z.read('xl/workbook.xml'))
    sheets=doc.findall('{*}sheets/{*}sheet')
    assert len(sheets)==len(list((source/'tables').glob('*.csv')))+1
    assert len({s.attrib['name'] for s in sheets})==len(sheets)

mi=rows(source/'tables/Supplementary_Clinical_C2_multiple_imputation_EN.csv')
assert len(mi)==9 and all(int(r['valid_imputations'])==50 for r in mi)
diagnostics=rows(source/'tables/Adult_clinical_MI_diagnostics_EN.csv')
assert all(float(r['Rhat'])<1.1 for r in diagnostics)
composition=rows(source/'tables/Composition_permutation_diagnostics_EN.csv')
assert len(composition)==6
assert all(float(r['vegan_kernel_error'])<1e-7 and float(r['full_original_vegan_error'])<1e-7 for r in composition)
report={'status':'PASS','figures':len(figures),'figure_files':len(figures)*4,'workbook_sheets':len(sheets),'viewer_files_checked':len(rows(viewer/'file_manifest.csv')),'clinical_MI_models':len(mi),'composition_numerical_checks':len(composition),'raster_geometry':geometry}
(source/'audit/delivery_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='raster_geometry'},indent=2))
