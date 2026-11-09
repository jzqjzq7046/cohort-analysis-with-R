"""Read-only PPT inspection; copied media and text are technical audit artifacts."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import posixpath
import re
import xml.etree.ElementTree as ET

ROOT = Path('E:/YFS/workspace/FattyLiver_generation_analysis')
PPT = ROOT / '11.4 fatty liver .pptx'
LECTURE = Path('C:/Users/jzqjzq/Desktop/YFS/第六次课-医学统计报告规范解读.pptx')
OUT = ROOT / 'outputs/_directory_audit/statistical_reporting_review_20261105'
OUT.mkdir(parents=True, exist_ok=True)
ns = {'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
      'p': 'http://schemas.openxmlformats.org/presentationml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}

def md5(p):
    with open(p, 'rb') as f:
        return hashlib.file_digest(f, 'md5').hexdigest()

before = {str(p): md5(p) for p in (PPT, LECTURE)}
slides = []
with ZipFile(PPT) as z:
    names = sorted((n for n in z.namelist() if re.fullmatch(r'ppt/slides/slide\d+\.xml', n)),
                   key=lambda n: int(re.search(r'slide(\d+)', n).group(1)))
    for name in names:
        number = int(re.search(r'slide(\d+)', name).group(1))
        tree = ET.fromstring(z.read(name))
        texts = [v.text or '' for v in tree.findall('.//a:t', ns)]
        rel = 'ppt/slides/_rels/' + name.rsplit('/', 1)[1] + '.rels'
        rels = {v.attrib['Id']: v.attrib.get('Target', '') for v in ET.fromstring(z.read(rel))}
        images = []
        for j, pic in enumerate(tree.findall('.//p:pic', ns), 1):
            blip = pic.find('.//a:blip', ns)
            rid = blip.attrib.get('{' + ns['r'] + '}embed')
            target = posixpath.normpath(posixpath.join('ppt/slides', rels[rid]))
            data = z.read(target)
            path = OUT / f'slide_{number:02d}_image_{j}{Path(target).suffix}'
            path.write_bytes(data)
            images.append({'path': str(path), 'source_entry': target,
                           'md5': hashlib.md5(data).hexdigest(), 'bytes': len(data)})
        slides.append({'slide': number, 'text': '\n'.join(texts), 'images': images})
after = {str(p): md5(p) for p in (PPT, LECTURE)}
assert before == after
manifest = {'source_hashes_before': before, 'source_hashes_after': after,
            'inputs_unchanged': True, 'slides': slides}
(OUT / 'ppt_inspection.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
(OUT / 'ppt_text.txt').write_text('\n\n'.join(f"SLIDE {v['slide']}\n{v['text']}" for v in slides), encoding='utf-8')
for v in slides:
    print(json.dumps({'slide': v['slide'], 'images': v['images']}, ensure_ascii=False))
print('INPUTS_UNCHANGED', before)
