"""Verify completed delivery and save a durable completion record. Never shuts down."""
from pathlib import Path
from datetime import datetime, timezone
import json,csv,hashlib,shutil,math,xml.etree.ElementTree as ET
from PIL import Image
from pypdf import PdfReader
P=Path(__file__).resolve().parents[1]
def current(module):
 runs=[r for r in (P/'outputs'/module/'current').iterdir() if r.is_dir()]
 if len(runs)!=1 or (runs[0]/'run_status.txt').read_text(encoding='utf-8-sig').splitlines()[0]!='COMPLETED':raise ValueError(module)
 return runs[0]
def readcsv(path):
 with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def md5(path):return hashlib.md5(path.read_bytes()).hexdigest()
def main():
 modules=['06_Cross_generation_clinical','07_Cross_generation_microbiome_similarity','08_Microbiome_steatosis_diversity','09_Technical_covariates','10_Results_visualization','11_Nature_composites']
 runs=[current(m) for m in modules]
 for r in runs:
  if not any(v in r.name for v in ('approved_design','approved_results_summary')):raise ValueError('Superseded run: '+str(r))
 validation=readcsv(P/'approved_result_validation.csv')
 if not validation or any(r['passed'].upper()!='TRUE' for r in validation):raise ValueError('Statistical validation failed')
 figure=runs[-1];registry=readcsv(figure/'reports/figure_registry.csv')
 if len(registry)!=9:raise ValueError('Expected 5 main and 4 Supplementary figures')
 audit_rows=[]
 for r in registry:
  name=r['figure'];audit=figure/'audit';pdf=figure/'figures'/(name+'.pdf')
  for suffix in ('.alignment.json','.effective-text-audit.json','.collision-audit.json'):
   report=json.loads((audit/(name+suffix)).read_text(encoding='utf-8'))
   # The alignment engine records summary.verdict instead of top-level verdict in some versions.
   verdict=report.get('verdict',report.get('summary',{}).get('verdict'))
   if verdict!='PASS':raise ValueError(name+suffix+': '+str(verdict))
   if suffix=='.effective-text-audit.json' and report['minimum_pt']<4.9999:raise ValueError('Small text')
   if suffix=='.collision-audit.json' and report.get('summary',{}).get('fail',0):raise ValueError('Collision failures')
  page=PdfReader(pdf).pages[0]
  width=float(page.mediabox.width)*25.4/72;height=float(page.mediabox.height)*25.4/72
  # R's Cairo device truncates its PDF MediaBox to integer PostScript points.
  # Verify that exact device result, rather than expecting fractional-point
  # dimensions unsupported by this backend. The canvas remains within 183 mm.
  expected_width_pt=math.floor(float(r['nominal_width_mm'])*72/25.4)
  expected_height_pt=math.floor(float(r['nominal_height_mm'])*72/25.4)
  if (abs(float(page.mediabox.width)-expected_width_pt)>.001 or
      abs(float(page.mediabox.height)-expected_height_pt)>.001 or
      width>183.001 or height>170.001):raise ValueError('Nature/Cairo dimensions')
  svg=ET.parse(figure/'figures'/(name+'.svg'))
  if not any(e.tag.rsplit('}',1)[-1]=='text' for e in svg.iter()):raise ValueError('No editable SVG text')
  with Image.open(figure/'figures'/(name+'.png')) as image:
   dpi=image.info.get('dpi',(0,0))
   if min(dpi)<599:raise ValueError('Raster resolution')
  audit_rows.append(dict(figure=name,width_mm=width,height_mm=height,dpi=dpi,editable_svg=True,qa='PASS'))
 # Verify that all previous canonical figure files survive publication unchanged.
 old=json.loads((P/'preserved_previous_figure_manifest.json').read_text(encoding='utf-8'))
 for row in old:
  options=[P/'outputs'/row['module']/kind/row['run']/row['relative'] for kind in ('current','history')]
  found=[f for f in options if f.is_file()]
  if len(found)!=1 or md5(found[0])!=row['md5']:raise ValueError('Prior figure changed/lost: '+str(row))
 # Preserve the first-pass auxiliary drafts outside the public scientific figure gallery.
 # Statistical CSVs/RDSs are untouched; final manuscript figures are the 11-module exports.
 div=runs[2];drafts=div/'audit/initial_style_drafts';drafts.mkdir(exist_ok=True)
 moved=[]
 for stem in ('Figure2_Shannon_distribution','Figure3_Shannon_adjusted_difference'):
  for ext in ('.pdf','.png','.svg'):
   src=div/'diversity'/(stem+ext)
   if src.exists():
    dst=drafts/src.name
    if dst.exists():raise ValueError('Draft archive collision')
    if not src.resolve().is_relative_to(P.resolve()) or not dst.resolve().is_relative_to(P.resolve()):raise ValueError('Unsafe draft path')
    checksum=md5(src);shutil.move(str(src),str(dst))
    if md5(dst)!=checksum:raise ValueError('Draft archive copy changed')
    moved.append(dict(file=src.name,md5=checksum))
 (drafts/'README_zh.md').write_text('这些是首次统计进程中已加载旧绘图函数生成的辅助草图。统计结果未改变；正式配图由11独立入口按新设计重绘并完成QA。为避免旧模板图注混入正式查看目录，将草图原样保留于audit。\n',encoding='utf-8')
 (drafts/'preserved_files.json').write_text(json.dumps(moved,indent=2),encoding='utf-8')
 summary=runs[4]/'reports/results_summary_zh.md';shutil.copy2(summary,P/'outputs/RESULTS_zh.md')
 for r in (runs[4],runs[5]):
  shutil.copy2(P/'approved_result_validation.csv',r/'audit/approved_result_validation.csv')
  shutil.copy2(P/'preserved_previous_figure_manifest.json',r/'audit/preserved_previous_figure_manifest.json')
  files=['config/approved_design_20261020.md','R/common.R','R/approved_design.R','R/approved_visualization.R','R/nature_export_approved.R']+[f'scripts/{s}' for s in ('06_cross_generation_clinical.R','07_cross_generation_microbiome_similarity.R','08_microbiome_steatosis_diversity.R','09_joint_technical_covariates.R','10_visualize_results.R','11_nature_composites.R')]
  with (r/'audit/final_code_fingerprints.csv').open('w',encoding='utf-8-sig',newline='') as f:
   w=csv.DictWriter(f,fieldnames=['file','md5']);w.writeheader();w.writerows(dict(file=name,md5=md5(P/name)) for name in files)
 readme=P/'README.md'
 if readme.exists():
  text=readme.read_text(encoding='utf-8-sig');marker='## 确认设计重分析完成（2026-10-20）'
  if marker not in text:
   text=marker+'\n\n当前模型与展示以[确认设计](config/approved_design_20261020.md)为准；[中文结果摘要](outputs/RESULTS_zh.md)、[重点表格](Results_表格与图片/Results_重点表格/index.html)、[Nature主要及补充合图](outputs/11_Nature_composites/index.html)。旧结果/代码保持可恢复。模型不再默认增加药物或营养层，未成年人保留。运行说明：[中文](README_approved_design_zh.md) / [English](README_approved_design_en.md)。下方早期更新记录为历史。\n\n'+text
   readme.write_text(text,encoding='utf-8')
 result=dict(status='COMPLETED',client_date='2026-10-20',written_at_utc=datetime.now(timezone.utc).isoformat(),
  design='approved_design_v1_20261020',shutdown_authorized=True,validation_passed=True,validation_checks=len(validation),figure_qa_passed=True,
  main_figures=5,supplementary_figures=4,figure_files=27,prior_figures_preserved=len(old),completed_runs=[str(r) for r in runs],
  figures=audit_rows,summary=str(summary),key_tables=str(P/'Results_表格与图片/Results_重点表格/index.html'),
  shutdown_note='Authorized by user; separate helper will request non-forced shutdown only after this record is verified.')
 (P/'completion_approved_design_20261020.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
 print('DELIVERY VERIFIED:',len(validation),'statistical checks;',len(registry),'figures;',len(old),'prior figure files preserved')
 print('Completion record:',P/'completion_approved_design_20261020.json')
if __name__=='__main__':main()
