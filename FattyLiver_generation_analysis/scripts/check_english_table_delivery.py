"""Final file/link/provenance checks; no statistical analyses."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
import csv, hashlib, json
from consolidate_key_table_folder import rebuild_index

ROOT=Path(__file__).resolve().parents[1]
class Links(HTMLParser):
 def __init__(self):super().__init__();self.links=[]
 def handle_starttag(self,tag,attrs):
  if tag=='a':self.links += [v for k,v in attrs if k=='href']

def main():
 rebuild_index()
 runs=[p for p in (ROOT/'outputs/15_English_journal_tables/current').iterdir() if p.is_dir()]
 assert len(runs)==1
 run=runs[0];checks=[]
 def check(name,ok):
  assert ok,name
  checks.append(dict(check=name,passed=True))
 check('Legacy root directory deleted',not (ROOT/'Results_重点表格').exists())
 base=ROOT/'Results_表格与图片/Results_重点表格';view=base/'tables/English'
 check('Canonical key-table folder exists',base.is_dir())
 with (view/'file_manifest.csv').open(encoding='utf-8-sig',newline='') as f:
  for r in csv.DictReader(f):
   a=view/r['file'];b=Path(r['source'])
   check('Published copy matches '+r['file'],hashlib.md5(a.read_bytes()).hexdigest()==r['md5']==hashlib.md5(b.read_bytes()).hexdigest())
 for p in [base/'index.html']+list(view.glob('*.html')):
  parser=Links();parser.feed(p.read_text(encoding='utf-8'))
  for link in parser.links:
   parsed=urlsplit(link)
   if parsed.scheme or not parsed.path:continue
   check('Local link '+p.name+' -> '+link,(p.parent/unquote(parsed.path)).exists())
 qa=ROOT/'inventory/key_table_folder_migration_20261020/final_predeletion_checks.json'
 for r in json.loads(qa.read_text(encoding='utf-8-sig')):
  check('Preserved legacy file '+r['relative'],hashlib.sha256(Path(r['preserved_at']).read_bytes()).hexdigest().upper()==r['sha256'])
 with (run/'audit/validation_checks.csv').open(encoding='utf-8-sig',newline='') as f:
  old=list(csv.DictReader(f));check('All 930 workbook/source checks passed',len(old)==930 and all(r['passed']=='True' for r in old))
 record=dict(canonical=str(base),english_tables=str(view),original_root_deleted=True,checks=checks,
  visual_review='All ten sheets reviewed. Final changes to numeric alignment and P precision reviewed in Clinical, Shannon and Shannon sensitivity previews. No clipping, error cells or missing-as-zero substitution identified.',
  scope='Presentation-only export. Original analytical results, raw data and figures preserved. Original adult complete-case Table S1 explicitly separated from later MI.')
 (run/'audit/final_delivery_checks.json').write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
 (ROOT/'completion_english_tables_20261020.json').write_text(json.dumps({k:v for k,v in record.items() if k!='checks'},ensure_ascii=False,indent=2),encoding='utf-8')
 print('PASS:',len(checks),'delivery, link and migration checks; 930 workbook/source checks passed')
 print('English tables:',view)

if __name__=='__main__':main()
