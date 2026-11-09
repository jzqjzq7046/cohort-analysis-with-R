"""Publish approved presentation edits; no statistical analysis is executed."""
from pathlib import Path
import csv,hashlib,html,json,shutil,subprocess
from openpyxl import load_workbook  # read-only verification; authoring uses Artifact Tool
root=Path(__file__).resolve().parents[1]
run=(root/'outputs/24_Two_generation_reanalysis/LATEST_RUN.txt').read_text().strip()
current=root/'outputs/24_Two_generation_reanalysis/current'/run
source=current/'tables/Manuscript_Table1_Table2_20261024';audit=current/'audit/Table1_Table2_20261024'
spec=json.loads((audit/'table_spec.json').read_text(encoding='utf-8'))['Table2']
cols=spec['columns'];assert len(cols)==8 and len(spec['rows'])==6
assert cols[2]=='Participants of younger generation, N'
assert {r[1] for r in spec['rows']}=={'Unadjusted','Model 1'}
raw=list(csv.DictReader((current/'tables/Table2_Clinical_cross_generation_all_models_EN.csv').open(encoding='utf-8-sig')))
for row in spec['rows']:
 match=[r for r in raw if r['comparison']==row[0] and r['model']==('C1' if row[1]=='Model 1' else 'Unadjusted')]
 assert len(match)==1
 r=match[0]
 assert int(r['n'])==row[2] and abs(float(r['OR'])-row[5])<1e-12 and abs(float(r['p_value'])-row[7])<1e-14
 assert row[6]==f"{float(r['CI_low']):.2f}–{float(r['CI_high']):.2f}"
wb=load_workbook(source/'Table1_Table2_EN.xlsx',data_only=True)
assert wb.sheetnames==['Table1','Table2'] and wb['Table1']['B5'].value==5660
ws=wb['Table2'];assert ws['D4'].value==spec['group_header']
for i,row in enumerate(spec['rows'],6):
 for j,value in enumerate(row,1):
  actual=ws.cell(i,j).value
  assert abs(actual-value)<1e-11 if isinstance(value,(int,float)) else actual==value
assert 'D4:E4' in {str(m) for m in ws.merged_cells.ranges}
css='body{font:15px/1.55 "Times New Roman",serif;margin:30px;color:#111}table{border-collapse:collapse;width:100%;border-top:1.5px solid #111;border-bottom:1.5px solid #111}th{padding:8px;text-align:center}thead{border-bottom:1px solid #111}thead tr:first-child th[colspan]{border-bottom:1px solid #777}td{padding:6px 8px;text-align:center}.notes{font-size:13px;margin-top:14px}a{color:#245882}.scroll{overflow:auto}h2{text-align:center;font-size:20px}@media print{body{margin:0}nav{display:none}}'
lines=['<!doctype html><meta charset="utf-8">',f'<style>{css}</style>',f'<h2>{html.escape(spec["title"])}</h2>','<div class="scroll"><table><thead><tr>']
for j,col in enumerate(cols):
 if j==3:lines.append(f'<th colspan="2">{html.escape(spec["group_header"])}</th>')
 elif j!=4:lines.append(f'<th rowspan="2">{html.escape(col)}</th>')
lines.append('</tr><tr>'+''.join(f'<th>{html.escape(cols[j])}</th>' for j in [3,4])+'</tr></thead><tbody>')
for row in spec['rows']:
 values=[]
 for j,v in enumerate(row):
  if isinstance(v,(int,float)):
   if j==2:v=f'{v:,.0f}'
   elif j==7:v=f'{v:.2e}' if v<.0001 else f'{v:.4f}'
   else:v=f'{v:.2f}'
  values.append(html.escape(str(v)))
 lines.append('<tr>'+''.join(f'<td>{v}</td>' for v in values)+'</tr>')
lines+=['</tbody></table></div><div class="notes">']+['<p>'+html.escape(n)+'</p>' for n in spec['notes']]+['</div>','<nav><a href="Table2_Clinical_association.pdf">PDF</a> | <a href="Table2_Clinical_association.png">PNG</a> | <a href="Table2_Clinical_association.tiff">TIFF</a> | <a href="Table2_Clinical_association.svg">SVG</a> | <a href="Table2_Clinical_association.csv">CSV</a> | <a href="Table1_Table2_EN.xlsx">Excel</a></nav>']
(source/'Table2_Clinical_association.html').write_text('\n'.join(lines),encoding='utf-8')
for support in source.glob('*.ndjson'):shutil.move(str(support),str(audit/support.name))
report={'status':'PASS','scope':'Table 2 presentation only','rows':6,'all_estimates_match_saved_models':True,'Table1_N':5660,'approved_group_header':spec['group_header'],'Excel_merged_header':'D4:E4','models':['Unadjusted','Model 1']}
(audit/'presentation_update_20261030.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(report,indent=2,ensure_ascii=False))
subprocess.run(['C:/Users/jzqjzq/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe',str(root/'tools/organize_results.py'),'--root',str(root)],check=True)
