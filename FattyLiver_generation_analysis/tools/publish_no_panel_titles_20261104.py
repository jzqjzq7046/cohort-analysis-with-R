"""Publish only the five reviewed figure sets and replace their PPT media bytes."""
from pathlib import Path
import csv, hashlib, json, posixpath, re, shutil, sys, zipfile
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
AUDIT=ROOT/'inventory/no_panel_titles_20261104'
TARGETS=json.loads((AUDIT/'targets.json').read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def stems(t):return t['stems'] if isinstance(t['stems'],list) else [t['stems']]
def check_statistics():
    with (AUDIT/'unchanged_statistics.csv').open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    assert rows
    for r in rows:
        assert hashlib.md5(Path(r['file']).read_bytes()).hexdigest()==r['md5'],r['file']
    return len(rows)
def figures():
    sys.path.insert(0,str(ROOT/'.figure_qa_python'))
    import pymupdf as fitz
    from PIL import Image
    before={str(f):sha(f) for t in TARGETS for f in (Path(t['parent'])/'figures').rglob('*') if f.is_file()}
    allowed=[]
    removed=['Observed prevalence','Cross-generation association','Observed similarity',
             'Conditional similarity','Observed diversity','B1 association','Model 3 association',
             'Overall steatosis associations','Overall steatosis-associated abundance']
    for t in TARGETS:
        src=AUDIT/'staging'/t['module']; dest=Path(t['parent'])
        for stem in stems(t):
            with fitz.open(src/'figures'/f'{stem}.pdf') as doc:
                text=' '.join(page.get_text() for page in doc)
            assert not any(s in text for s in removed),stem
            for ext in ('pdf','png','tiff','svg'):
                f=src/'figures'/f'{stem}.{ext}'
                assert f.exists() and f.stat().st_size>1000
                if ext in ('png','tiff'):
                    with Image.open(f) as im:assert all(abs(float(v)-600)<.1 for v in im.info['dpi'])
                target=dest/'figures'/f.name
                shutil.copyfile(f,target);assert sha(f)==sha(target);allowed.append(str(target))
            for f in (src/'audit').glob(stem+'.*'):
                target=dest/'audit'/f.name
                backup=AUDIT/'previous_figures'/t['module']/'audit'/f.name
                if target.exists() and not backup.exists():
                    backup.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(target,backup)
                shutil.copyfile(f,target)
    changed=[f for f,h in before.items() if sha(f)!=h]
    assert set(changed)==set(allowed)
    checked=check_statistics()
    report=dict(status='PASS',figures=5,exports=20,only_requested_figure_files_changed=True,
                protected_statistical_files_unchanged=checked,changed=changed)
    (AUDIT/'published_figures.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='changed'}))
def pptx():
    source=ROOT/'11.4 fatty liver .pptx'
    frozen_hash=sha(source)
    backup=AUDIT/'11.4_before_five_title_removals.pptx'
    if backup.exists():assert sha(backup)==frozen_hash,'User source changed since the saved backup'
    else:shutil.copy2(source,backup)
    mapping={s:Path(t['parent'])/'figures'/f'{s}.tiff' for t in TARGETS for s in stems(t)}
    ns={'p':'http://schemas.openxmlformats.org/presentationml/2006/main',
        'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
        'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
    changes={};found=[];all_images=[]
    with zipfile.ZipFile(source) as z:
        for n in z.namelist():
            if not re.fullmatch(r'ppt/slides/slide\d+[.]xml',n):continue
            tree=ET.fromstring(z.read(n));text=' '.join(tree.itertext())
            relpath=posixpath.dirname(n)+'/_rels/'+posixpath.basename(n)+'.rels'
            rels={r.get('Id'):r.get('Target') for r in ET.fromstring(z.read(relpath))}
            pics=tree.findall('.//p:pic',ns)
            for pic in pics:
                desc=pic.find('p:nvPicPr/p:cNvPr',ns).get('descr','')
                rid=pic.find('.//a:blip',ns).get('{'+ns['r']+'}embed')
                media=posixpath.normpath(posixpath.join(posixpath.dirname(n),rels[rid]));all_images.append(media)
                stem=desc if desc in mapping else None
                if not stem and not desc and len(pics)==1 and '16 SGBs' in text and 'abundance' in text:
                    stem='Figure25_Overall_abundance_ratios'
                if stem:
                    ext=Path(media).suffix.lower()
                    assert ext in ('.png','.tif','.tiff'),media
                    asset=mapping[stem].with_suffix('.tiff' if ext in ('.tif','.tiff') else ext)
                    changes[media]=asset.read_bytes()
                    found.append(dict(slide=n,figure=stem,media=media))
        assert {r['figure'] for r in found}==set(mapping) and len(found)==5,found
        assert all(all_images.count(m)==1 for m in changes),'Shared media requires separate relationship'
        candidate=AUDIT/'11.4_five_titles_removed.pptx'
        with zipfile.ZipFile(candidate,'w') as out:
            for info in z.infolist():out.writestr(info,changes.get(info.filename,z.read(info.filename)))
    with zipfile.ZipFile(backup) as old,zipfile.ZipFile(candidate) as new:
        changed=[n for n in old.namelist() if old.read(n)!=new.read(n)]
        assert set(changed)==set(changes) and len(changed)==5
        assert new.testzip() is None
    assert sha(source)==frozen_hash,'User changed the presentation during this edit'
    shutil.copyfile(candidate,source);assert sha(source)==sha(candidate)
    report=dict(status='PASS',source_sha256=frozen_hash,updated_sha256=sha(source),
        replaced_images=found,changed_parts=changed,all_slide_xml_and_captions_unchanged=True)
    (AUDIT/'ppt_media_replacement.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
{'figures':figures,'pptx':pptx}[sys.argv[1]]()
