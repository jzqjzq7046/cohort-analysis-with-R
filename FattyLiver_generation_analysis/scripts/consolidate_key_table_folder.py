"""Merge the legacy viewing folder into the user-designated canonical folder.

Copy and checksum first. This script never deletes either directory.
"""
from pathlib import Path
from datetime import datetime
import json, hashlib, shutil, html

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'Results_重点表格'
NEW=ROOT/'Results_表格与图片/Results_重点表格'
AUDIT=ROOT/'inventory/key_table_folder_migration_20261020'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 root=ROOT.resolve();old=OLD.resolve();new=NEW.resolve()
 if old!=root/'Results_重点表格' or new!=root/'Results_表格与图片/Results_重点表格':raise ValueError('Unexpected migration path')
 if not old.is_relative_to(root) or not new.is_relative_to(root) or old==new:raise ValueError('Unsafe migration')
 AUDIT.mkdir(parents=True,exist_ok=True);NEW.mkdir(parents=True,exist_ok=True)
 if not OLD.exists():
  print('Legacy folder already absent; canonical folder retained');return
 stamp=datetime.now().strftime('%Y%m%d_%H%M%S');rows=[]
 # Keep every unique file from both roots; colliding content is preserved in audit.
 for src in sorted(OLD.rglob('*')):
  if not src.is_file() or src.name.startswith('~$'):continue
  if src.is_symlink() or not src.resolve().is_relative_to(old):raise ValueError('Unexpected link in legacy directory')
  rel=src.relative_to(OLD);dst=NEW/rel;dst.parent.mkdir(parents=True,exist_ok=True)
  if dst.exists() and digest(src)!=digest(dst):
   backup=AUDIT/'premerge_destination'/rel;backup.parent.mkdir(parents=True,exist_ok=True)
   if not backup.exists():shutil.copy2(dst,backup)
  if not dst.exists() or digest(src)!=digest(dst):shutil.copy2(src,dst)
  if digest(src)!=digest(dst):raise AssertionError('Migration copy mismatch')
  rows.append(dict(relative=str(rel),sha256=digest(src),bytes=src.stat().st_size))
 (AUDIT/'copied_legacy_files.json').write_text(json.dumps(dict(old=str(old),canonical=str(new),files=rows),ensure_ascii=False,indent=2),encoding='utf-8')
 # Fix only browsing links that acquired one extra directory level.
 for p in NEW.rglob('*'):
  if p.is_file() and p.suffix.lower() in ('.html','.md'):
   text=p.read_text(encoding='utf-8-sig')
   updated=text.replace('../../outputs/','../../../outputs/').replace('../../config/','../../../config/')
   if updated!=text:p.write_text(updated,encoding='utf-8')
 # Index is a curated entry; copied historical indexes are preserved as files.
 index=NEW/'index.html'
 if index.exists():shutil.copy2(index,AUDIT/'legacy_root_index.html')
 rebuild_index()
 # Re-verify each copied file; revised relative-link documents have an audit copy.
 for r in rows:
  src=OLD/r['relative'];dst=NEW/r['relative']
  if digest(src)!=digest(dst):
   backup=AUDIT/'pre_link_revision'/r['relative'];backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,backup)
 (AUDIT/'copy_verification_complete.txt').write_text(f'{len(rows)} legacy files preserved in canonical folder or migration audit.\n',encoding='utf-8')
 patch_current_paths()
 print(f'Merged {len(rows)} legacy files; deletion may proceed after verification: {old}')

def rebuild_index():
 links=[]
 preferred=['tables/English/index.html','tables/index.html']
 for rel,label in zip(preferred,['English journal-style tables (recommended)','Original figure-to-table exports']):
  if (NEW/rel).exists():links.append(f'<li><a href="{rel}">{label}</a></li>')
 for p in sorted(NEW.glob('adult_MI_*/index.html'),reverse=True):links.append(f'<li><a href="{p.parent.name}/index.html">Adult multiple imputation — {p.parent.name}</a></li>')
 for p in sorted(NEW.glob('approved_design_*/index.html'),reverse=True):links.append(f'<li><a href="{p.parent.name}/index.html">Approved primary results — {p.parent.name}</a></li>')
 histories=sorted([p for p in NEW.glob('20*/index.html')],reverse=True)
 history=''.join(f'<li><a href="{p.parent.name}/index.html">{p.parent.name}</a></li>' for p in histories)
 page='<!doctype html><html lang="en"><meta charset="utf-8"><title>Key tables</title><style>body{font:16px/1.7 Arial;max-width:1000px;margin:36px auto;color:#222}a{color:#235584}</style><h1>Key tables</h1><p>Canonical viewing folder. Adopted analysis outputs remain in the project outputs directory.</p><ul>'+''.join(links)+'</ul><details><summary>Earlier result snapshots</summary><ul>'+history+'</ul></details></html>'
 (NEW/'index.html').write_text(page,encoding='utf-8')
 if (NEW/'tables/English/index.html').exists():
  (NEW/'CURRENT.txt').write_text('tables/English/index.html\n',encoding='utf-8')
  (NEW/'README_EN.md').write_text('# Canonical key-table viewer\n\nOpen index.html for the current English journal tables, original figure-to-table exports, adult multiple-imputation results, and earlier snapshots.\n\nOnly this nested key-table folder is used for future exports. Authoritative analysis runs remain in outputs; this folder contains viewing copies. The English Table S1 is explicitly the original adult complete-case version; the later MI results are a separate entry.\n',encoding='utf-8')

def patch_current_paths():
 files=list((ROOT/'scripts').glob('*.py'))+list((ROOT/'scripts').glob('*.R'))+list((ROOT/'scripts').glob('*.md'))+[ROOT/'README.md']
 for p in files:
  if p.name==Path(__file__).name:continue
  text=p.read_text(encoding='utf-8-sig');updated=text
  updated=updated.replace("ROOT/'Results_重点表格", "ROOT/'Results_表格与图片/Results_重点表格")
  updated=updated.replace("P/'Results_重点表格", "P/'Results_表格与图片/Results_重点表格")
  updated=updated.replace('file.path(root,"Results_重点表格/', 'file.path(root,"Results_表格与图片/Results_重点表格/')
  updated=updated.replace('](Results_重点表格/', '](Results_表格与图片/Results_重点表格/')
  updated=updated.replace('`Results_重点表格/tables`','`Results_表格与图片/Results_重点表格/tables`')
  # Only the changed operational files are backed up; historical code is untouched.
  if updated!=text:
   backup=AUDIT/'previous_code'/p.relative_to(ROOT);backup.parent.mkdir(parents=True,exist_ok=True)
   if not backup.exists():shutil.copy2(p,backup)
   p.write_text(updated,encoding='utf-8')
 # Gallery rebuilds preserve this nested folder and should expose its entry.
 p=ROOT/'R/result_gallery.R';text=p.read_text(encoding='utf-8-sig')
 needle="<nav><a href='#tables'>表格</a>"
 replacement="<nav><a href='Results_重点表格/index.html'>重点表格 / Key tables</a><a href='#tables'>表格</a>"
 if needle in text:
  backup=AUDIT/'previous_code/R/result_gallery.R';backup.parent.mkdir(parents=True,exist_ok=True)
  if not backup.exists():shutil.copy2(p,backup)
  p.write_text(text.replace(needle,replacement),encoding='utf-8')
 # Record durable instruction for future exporters.
 p=ROOT/'AGENTS.md';text=p.read_text(encoding='utf-8')
 rule='- The only key-table viewing folder is Results_表格与图片/Results_重点表格; never recreate Results_重点表格 at project root. English manuscript tables go in its tables/English subfolder.\n'
 if rule not in text:p.write_text(text+'\n'+rule,encoding='utf-8')

if __name__=='__main__':main()
