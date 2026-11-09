"""Verify and publish model-label-only redraws; no numerical fits or captions."""
from pathlib import Path
import csv,collections,hashlib,json,re,shutil,sys
root=Path(__file__).resolve().parents[1];audit=root/'inventory/model_labels_20261030'
sys.path.insert(0,str(root/'.figure_qa_python'));import pymupdf
from PIL import Image
assert (audit/'STATUS').read_text().strip()=='STAGED_RENDERING_COMPLETED'
protected=json.loads((audit/'protected_before.json').read_text())
def md5(p):return hashlib.md5(Path(p).read_bytes()).hexdigest()
assert all(md5(r['path'])==r['md5'] for r in protected),'A protected statistic or legend changed'
rows=list(csv.DictReader((root/'Results_表格与图片/file_manifest.csv').open(encoding='utf-8-sig')))
by_source={r['source_path']:r for r in rows}
targets=json.loads((audit/'targets.json').read_text())
raw=re.compile(r'(?<![A-Za-z0-9])(?:C1|C2|B1|B2|P1|P2)(?![A-Za-z0-9])');models=re.compile(r'Model\s*[1-6]')
def content(p):
 with pymupdf.open(p) as d:return '\n'.join(page.get_text() for page in d),[(page.rect.width,page.rect.height) for page in d]
def numbers(text):
 text=models.sub('',raw.sub('',text)).replace('−','-').replace('–','-')
 return collections.Counter(re.findall(r'(?<![A-Za-z])[-+]?(?:\d*\.\d+|\d+)(?:e[-+]?\d+)?',text))
checks=[]
for t in targets:
 old=root/'Results_表格与图片'/by_source[t['source_pdf']]['viewing_path']
 new=Path(t['source_pdf']) if t['stem']=='Study_population_and_analysis_sets' else audit/'staging'/t['module']/'figures'/(t['stem']+'.pdf')
 oldtext,oldsize=content(old);newtext,newsize=content(new)
 assert not raw.search(newtext),(t['stem'],'old labels remain')
 assert oldsize==newsize,(t['stem'],'canvas changed')
 a,b=numbers(oldtext),numbers(newtext)
 assert a==b,(t['stem'],'numerical text changed',dict(a-b),dict(b-a))
 if t['stem']=='Study_population_and_analysis_sets':qa=Path(t['source_pdf']).parents[2]/'audit/Cohort_overview_two_generation_20261024'
 else:qa=new.parent.parent/'audit'
 for suffix in ['collision-audit.json','effective-text-audit.json']:
  q=json.loads((qa/(t['stem']+'.'+suffix)).read_text());assert q['verdict']=='PASS',(t['stem'],suffix,q)
 for ext in ['.pdf','.png','.tiff','.svg']:assert new.with_suffix(ext).exists()
 for ext in ['.png','.tiff']:
  with Image.open(new.with_suffix(ext)) as im:assert min(im.info.get('dpi',(0,0)))>599
 checks.append(dict(module=t['module'],figure=t['stem'],new_pdf=str(new),old_pdf=str(old),old_labels_removed=True,same_canvas=True,same_numeric_text=True,rendered_QA='PASS'))
(audit/'prepublication_checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
# All checks pass before the first generic current figure is replaced.
copied=[]
for t,c in zip(targets,checks):
 old=Path(c['old_pdf']);new=Path(c['new_pdf']);dest=Path(t['source_pdf'])
 base=root/'outputs'/t['module'];run=(base/'LATEST_RUN.txt').read_text().strip();cur=base/'current'/run
 backup=cur/'audit/model_labels_20261030/before';backup.mkdir(parents=True,exist_ok=True)
 for ext in ['.pdf','.png','.tiff','.svg']:
  src=new.with_suffix(ext);target=dest.with_suffix(ext);prior=old.with_suffix(ext)
  saved=backup/(t['stem']+ext)
  if not saved.exists():shutil.copy2(prior,saved)
  if src.resolve()!=target.resolve():shutil.copy2(src,target)
  assert md5(src)==md5(target)
  copied.append(str(target))
 if t['stem']!='Study_population_and_analysis_sets':
  qa_dest=cur/'audit/model_labels_20261030/final';qa_dest.mkdir(parents=True,exist_ok=True)
  for f in (new.parent.parent/'audit').glob(t['stem']+'.*'):shutil.copy2(f,qa_dest/f.name)
assert all(md5(r['path'])==r['md5'] for r in protected)
# Audit every current figure, including those that needed no redraw.
all_checked=0
for r in rows:
 if r['viewing_path'].startswith('01_当前分析/') and r['category']=='Figures_图片' and r['file'].endswith('.pdf'):
  text,_=content(r['source_path']);assert not raw.search(text),(r['file'],'stale model text');all_checked+=1
report={'status':'PASS','figures_redrawn':len(checks),'current_figure_PDFs_checked':all_checked,'four_format_files_updated':len(copied),'protected_inference_and_caption_files':len(protected),'inference_and_captions_unchanged':True,'mapping':{'C1':'Model 1','C2':'Model 2','B1':'Model 3','B2':'Model 4','P1':'Model 5','P2':'Model 6'},'unadjusted_clinical_label':'Unadjusted'}
(audit/'publication_validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
