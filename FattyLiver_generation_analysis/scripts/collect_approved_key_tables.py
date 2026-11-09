"""Create a versioned, checksum-verified key-table snapshot; never fit models."""
from pathlib import Path
from datetime import datetime
import csv, hashlib, html, shutil
ROOT=Path(__file__).resolve().parents[1]
def current(module, approved=False):
 runs=[p for p in (ROOT/'outputs'/module/'current').iterdir() if p.is_dir()]
 if len(runs)!=1 or (runs[0]/'run_status.txt').read_text(encoding='utf-8-sig').splitlines()[0]!='COMPLETED':raise ValueError(module)
 if approved and 'approved_design_v1' not in runs[0].name:raise ValueError('Superseded source: '+module)
 return runs[0]
def md5(p):return hashlib.md5(p.read_bytes()).hexdigest()
def main():
 dest=ROOT/'Results_表格与图片/Results_重点表格'/('approved_design_'+datetime.now().strftime('%Y%m%d_%H%M%S'));dest.mkdir(parents=True)
 clinical=current('06_Cross_generation_clinical',True)/'clinical'
 similarity=current('07_Cross_generation_microbiome_similarity',True)/'similarity'
 diversity=current('08_Microbiome_steatosis_diversity',True)/'diversity'
 sources=[(current('01_Table1')/'tables/Table1_by_generation.csv','01_Table1_三代基本情况.csv'),
  (current('10_Results_visualization')/'tables/Main_clinical_comparison_wide.csv','02_跨代脂肪肝_粗与主要调整.csv'),
  (similarity/'Main_similarity_adjusted.csv','03_家系菌群相似性_主要调整.csv'),
  (diversity/'Main_Shannon_adjusted.csv','04_脂肪肝_Shannon主要调整.csv'),
  (diversity/'Main_composition_adjusted.csv','05_脂肪肝_整体组成主要调整.csv'),
  (clinical/'Supplementary_clinical_adult.csv','S01_临床成人扩展.csv'),
  (similarity/'Supplementary_similarity.csv','S02_相似性敏感性.csv'),
  (diversity/'review_Shannon_models.csv','S03_Shannon全部方案.csv'),
  (diversity/'review_composition_association.csv','S04_整体组成全部方案.csv'),
  (clinical/'review_same_sample_audit.csv','A01_临床同样本核查.csv'),
  (diversity/'review_sample_counts.csv','A02_菌群样本与病例数.csv')]
 manifest=[]
 for src,name in sources:
  dst=dest/name;shutil.copy2(src,dst)
  if md5(src)!=md5(dst):raise ValueError('Copy checksum mismatch')
  manifest.append(dict(file=name,source=str(src),source_md5=md5(src),copied_md5=md5(dst)))
 with (dest/'file_manifest.csv').open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(manifest[0]));w.writeheader();w.writerows(manifest)
 summary=current('10_Results_visualization')/'reports/results_summary_zh.md';shutil.copy2(summary,dest/'结果摘要_zh.md')
 shutil.copy2(ROOT/'config/approved_design_20261020.md',dest/'确认设计_zh.md')
 body='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>脂肪肝重点表格</title><body><h1>确认设计：重点表格</h1><p>主表01–05；补充S01–S04；A为审计。全部原始P。来源和MD5见file_manifest.csv。</p><ul>'
 body+=''.join('<li><a href="'+html.escape(r['file'],quote=True)+'">'+html.escape(r['file'])+'</a></li>' for r in manifest)
 body+='</ul><a href="结果摘要_zh.md">中文结果摘要</a> · <a href="确认设计_zh.md">设计</a></body></html>'
 (dest/'index.html').write_text(body,encoding='utf-8')
 (dest.parent/'LATEST.txt').write_text(dest.name+'\n',encoding='utf-8')
 (dest.parent/'index.html').write_text('<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0;url='+dest.name+'/index.html"><a href="'+dest.name+'/index.html">最新重点表格</a>',encoding='utf-8')
 from consolidate_key_table_folder import rebuild_index
 rebuild_index()
 print('Key tables:',dest);print('Verified unchanged copies:',len(sources))
if __name__=='__main__':main()
