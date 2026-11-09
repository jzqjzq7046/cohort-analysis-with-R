"""Prepare/verify/publish tables corresponding to the adopted nine Nature figures.

Only aggregate source data are read. No model is refitted, no participant data exported.
"""
from pathlib import Path
import csv, hashlib, html, json, math, shutil, sys, zipfile
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
PAIR=['G0-G1','G1-G2','G0-G2']
SCENARIOS={
 'core_base':'技术基础参照（同主要样本）','core_adjusted':'全龄主要调整',
 'adult_core':'成人核心参照','adult_lifestyle':'成人＋生活方式','adult_BMI':'成人＋生活方式及BMI',
 'bowel_base':'排便资料完整参照','bowel_adjusted':'＋排便及乳酸菌制剂',
 'unperturbed_core':'排除近期抗生素或胃肠炎','DNA_base':'DNA资料完整参照','core_DNA':'＋DNA浓度'}
SOURCE_NAMES=['age','proportions','clinical_prevalence_same_sample','clinical','distance','pair','shannon','composition','samples','heterogeneity']
GENERAL='范围固定为基线319个三代主家系。所有P为原始P。表格对应原9张图，不混入之后的成人多重插补结果。'
CLINICAL='核心：年轻年龄/性别，年长观察人数/平均年龄/女性比例，调查地区/城乡；恒定项移除。家系GEE。'
MICRO='双方年龄（自然样条）、性别、地区、建库批次、测序池、矩阵总量、采便季节、冷冻间隔、近期抗生素；含同地区/批次/池指标。'
DIVERSITY='年龄（各代自然样条）、性别、代际、地区、建库批次、测序池、矩阵总量、季节、冷冻间隔、采便至超声时间差、近期抗生素。'
SENS='成人核心→生活方式→BMI使用共同完整样本；排便、DNA各有同样本参照；排除近期扰动会改变样本。技术基础参照已调整部分变量。'

def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def number(row,key):
 x=row.get(key,'')
 if x in ('','NA','NaN',None):return None
 n=float(x)
 if not math.isfinite(n):raise ValueError(f'Nonfinite {key}')
 return int(n) if n.is_integer() else n
def md5(p):return hashlib.md5(p.read_bytes()).hexdigest()
def current():
 runs=[r for r in (ROOT/'outputs/11_Nature_composites/current').iterdir() if r.is_dir()]
 if len(runs)!=1 or (runs[0]/'run_status.txt').read_text(encoding='utf-8-sig').splitlines()[0]!='COMPLETED':raise ValueError('Completed figure run required')
 return runs[0]
def writecsv(path,headers,rows):
 with path.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.writer(f);w.writerow(headers);w.writerows(rows)
def make(ident,sheet,title,figure,panels,recommendation,headers,rows,formats,widths,notes,sources):
 assert all(len(r)==len(headers) for r in rows)
 return dict(id=ident,sheet=sheet,title=title,figure=figure,panels=panels,recommendation=recommendation,headers=headers,rows=rows,formats=formats,widths=widths,notes=notes,sources=sources)

def prepare(out):
 src=current();sd=src/'source_data';d={k:read(sd/(k+'.csv')) for k in SOURCE_NAMES}
 tables=[]
 ages={r['generation']:r for r in d['age']}
 rows=[]
 for r in d['proportions']:
  # The source carries a zero-person bookkeeping row absent from the figure.
  if r['group']=='G2 age unavailable' and number(r,'n')==0:continue
  group=r['group'];a=ages.get(group,{})
  rows.append([group,number(r,'n'),number(r,'cases'),number(r,'proportion')]+[number(a,k) for k in ['median','q1','q3','minimum','maximum']])
 tables.append(make('01','01_三代基本情况','三代样本、年龄及脂肪肝分布','Fig01_Cohort','b–d','部分转表；保留a流程图',
  ['分组','人数','脂肪肝病例','比例','年龄中位数（岁）','年龄Q1（岁）','年龄Q3（岁）','最小年龄（岁）','最大年龄（岁）'],rows,
  ['text','count','count','percent1','one','one','one','one','one'],[135,85,110,95,155,125,125,145,145],
  ['原始描述，不调整协变量。年龄Q1–Q3是四分位范围，不是95%CI。','G2年龄亚组包含在G2及Overall中，不能将所有行相加。','年龄仅在源表提供的G0/G1/G2行展示；未提供亚组年龄汇总的单元格留空。'],['age.csv','proportions.csv']))
 prev={(r['comparison'],r['older']):r for r in d['clinical_prevalence_same_sample']}
 rows=[]
 for pair in PAIR:
  r=next(r for r in d['clinical'] if r['comparison']==pair and r['model']=='R1_review_primary')
  a,b=prev[(pair,'No older case')],prev[(pair,'Older case')]
  rows.append([pair,number(a,'n'),number(a,'cases'),number(a,'cases')/number(a,'n'),number(b,'n'),number(b,'cases'),number(b,'cases')/number(b,'n'),number(r,'OR'),number(r,'CI_low'),number(r,'CI_high'),number(r,'p_value')])
 tables.append(make('02','02_脂肪肝跨代关联','脂肪肝跨代关联：观察比例及主要调整','Fig02_Clinical','a–b','适合整体转表',
  ['比较','无年长病例组N','无年长病例组病例','无年长病例组比例','有年长病例组N','有年长病例组病例','有年长病例组比例','调整OR','95%CI下限','95%CI上限','原始P'],rows,
  ['text','count','count','percent1','count','count','percent1','three','three','three','p'],[100,155,110,100,155,110,100,110,120,120,100],
  [CLINICAL,'观察比例与主要OR使用相同完整样本。G0–G1为317人/88例；涉及G2的全龄模型为484人/17例。','G1–G2的3.71是全龄主要OR，不能替换成之后成人插补分析的2.04。'],['clinical_prevalence_same_sample.csv','clinical.csv']))
 rows=[]
 for pair in PAIR:
  r=next(r for r in d['distance'] if r['comparison']==pair and r['scenario']=='core_adjusted')
  rows.append([pair]+[number(r,k) for k in ['n_people','families','within_pairs','between_pairs','observed_within_mean','observed_between_mean','adjusted_distance_difference','p_raw_two_sided']])
 tables.append(make('03','03_家系菌群相似性','家系内外菌群相似性：原始均值及调整差值','Fig03_Family_microbiome','a–b','适合整体转表',
  ['比较','人数','家系数','家系内配对数','家系间配对数','家系内距离均值','家系间距离均值','调整IF−BF差值','原始P'],rows,
  ['text']+['count']*4+['four']*3+['p'],[100,85,85,140,140,160,160,165,100],
  [MICRO,'IF=家系内，BF=家系间；距离越小越相似。人数可跨比较重复，配对数不是独立样本数。','9999次节点残差置换，最小add-one P为0.0001。没有独立配对95%CI，表中不制造CI。'],['distance.csv','pair.csv']))
 rows=[]
 for g in ['G0','G1','G2']:
  r=next(r for r in d['shannon'] if r['generation']==g and r['scenario']=='core_adjusted')
  rows.append([g]+[number(r,k) for k in ['n','cases','estimate','CI_low','CI_high','p_raw']])
 hetero=next(r for r in d['heterogeneity'] if r['scenario']=='core_adjusted')
 tables.append(make('04','04_Shannon主要结果','脂肪肝与Shannon多样性：主要调整结果','Fig04_Shannon','b','仅b转表；保留a原始分布',
  ['代际','人数','脂肪肝病例','调整Shannon差值','95%CI下限','95%CI上限','原始P'],rows,
  ['text','count','count','three','three','three','p'],[100,90,120,185,125,125,110],
  [DIVERSITY,'差值为脂肪肝减非脂肪肝；家系GEE及稳健95%CI。主要样本1267人/232例。',f'代际效应异质性原始P={number(hetero,"p_raw"):.4f}；不能仅凭G2点估计更负认定代际效应不同。'],['shannon.csv','samples.csv','heterogeneity.csv']))
 def comp_rows(data):
  rows=[]
  for scenario in SCENARIOS:
   for scheme in ['family','family_generation']:
    selected=[r for r in data if r['scenario']==scenario and r['scheme']==scheme]
    if not selected:continue
    r=selected[0];rows.append([SCENARIOS[scenario],'家系内' if scheme=='family' else '家系内同代']+[number(r,k) for k in ['n','cases','families','pseudo_F','R2_total_corrected','p_raw']])
  return rows
 tables.append(make('05','05_整体菌群组成','脂肪肝与整体菌群组成：主要条件关联','Fig05_Composition','b','仅b转表；保留a PCoA散点图',
  ['调整方案','置换范围','人数','脂肪肝病例','家系数','伪F','解释度R²（比例）','原始P'],comp_rows([r for r in d['composition'] if r['scenario']=='core_adjusted']),
  ['text','text','count','count','count','three','percent4','p'],[245,145,85,120,85,110,180,110],
  [DIVERSITY,'partial dbRDA额外条件控制家系。两种置换范围使用同一拟合，解释度相同并非重复样本。','Excel/HTML中的解释度按百分比显示；CSV保存0–1比例。P为9999次置换的F上尾原始P。'],['composition.csv']))
 models={'R2_adult_core_same_sample':'成人核心参照','R3_adult_lifestyle':'＋生活方式','R4_adult_lifestyle_BMI':'＋生活方式及BMI'}
 rows=[]
 for pair in PAIR:
  for model,label in models.items():
   r=next(r for r in d['clinical'] if r['comparison']==pair and r['model']==model)
   status='可估（稀少事件提示）' if r['status'].startswith('ESTIMATED') else '不可估：完全或近完全分离'
   rows.append([pair,label]+[number(r,k) for k in ['n','cases','OR','CI_low','CI_high','p_value']]+[status])
 tables.append(make('06','06_成人临床扩展原图','成人临床扩展：原图完整资料版本','Supp01_Clinical_adult','a–c','优先整体转表，清楚呈现不可估状态',
  ['比较','成人模型','人数','脂肪肝病例','OR','95%CI下限','95%CI上限','原始P','模型状态'],rows,
  ['text','text','count','count','three','three','three','p','text'],[100,185,85,120,105,120,120,110,295],
  [CLINICAL,'生活方式增加饮酒、吸烟、饮食评分、活动和能量；BMI层进一步增加BMI。三层使用共同完整样本。','严格对应原Supp01：G0–G1为256人/72例，涉及G2的成人模型为179人/10例；空数值表示不可估，不是0。','此表不包含后续MI；最新MI成人结果在重点表格另一个入口，不能把两个人群或版本混写。'],['clinical.csv']))
 rows=[]
 for scenario in SCENARIOS:
  for pair in PAIR:
   r=next(r for r in d['distance'] if r['comparison']==pair and r['scenario']==scenario)
   rows.append([SCENARIOS[scenario],pair]+[number(r,k) for k in ['n_people','families','adjusted_distance_difference','p_raw_two_sided']])
 tables.append(make('07','07_相似性敏感性','家系菌群相似性：全部敏感性方案','Supp02_Similarity_sensitivity','a–c','优先整体转表，便于横向比较样本与P',
  ['调整方案','比较','人数','家系数','调整IF−BF差值','原始P'],rows,
  ['text','text','count','count','four','p'],[290,100,100,100,180,120],
  [SENS,'主要集采用配对双方协变量；不加入采便–超声时间差。IF−BF为负表示家系内更相似。','不提供将重复配对当独立样本计算的CI。所有方案9999次受限节点置换，原始双侧P。'],['distance.csv']))
 rows=[]
 for scenario in SCENARIOS:
  for g in ['G0','G1','G2']:
   r=next(r for r in d['shannon'] if r['generation']==g and r['scenario']==scenario)
   rows.append([SCENARIOS[scenario],g]+[number(r,k) for k in ['n','cases','estimate','CI_low','CI_high','p_raw']])
 tables.append(make('08','08_Shannon敏感性','脂肪肝与Shannon：全部敏感性方案','Supp03_Shannon_sensitivity','a–c','优先整体转表，完整保留CI及P',
  ['调整方案','代际','人数','脂肪肝病例','调整差值','95%CI下限','95%CI上限','原始P'],rows,
  ['text','text','count','count','three','three','three','p'],[290,85,85,120,135,125,125,110],
  [SENS,'差值为脂肪肝减非脂肪肝；Gaussian家系GEE、稳健逐项95%CI及原始双侧P。','不同样本之间的效应变化不能全部归因于增加协变量。这里的成人方案是原图完整资料分析。'],['shannon.csv','samples.csv']))
 tables.append(make('09','09_组成敏感性','脂肪肝与整体菌群组成：全部敏感性方案','Supp04_Composition_sensitivity','a–b','优先整体转表，解释度与两种P同时展示',
  ['调整方案','置换范围','人数','脂肪肝病例','家系数','伪F','解释度R²（比例）','原始P'],comp_rows(d['composition']),
  ['text','text','count','count','count','three','percent4','p'],[290,145,85,120,85,100,180,110],
  [SENS,'同方案两种置换范围使用同一拟合，解释度相同但P可不同。CSV中的R²为0–1比例。','partial dbRDA控制家系及相应协变量；P为F上尾原始置换P，不把统计显著解释为解释度很大。'],['composition.csv']))
 assert [len(t['rows']) for t in tables]==[6,3,3,3,2,9,30,30,20]
 protected=list((src/'figures').glob('*'))+[sd/(k+'.csv') for k in SOURCE_NAMES]+[src/'reports/figure_legends.md']
 hashes=[dict(path=str(p),md5=md5(p)) for p in protected if p.is_file()]
 (out/'audit/source_fingerprints.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2),encoding='utf-8')
 for t in tables:
  t['csv']=t['id']+'_'+t['sheet'].split('_',1)[1]+'.csv'
  writecsv(out/'tables'/t['csv'],t['headers'],t['rows'])
 payload=dict(source_run=src.name,general=GENERAL,tables=tables,workbook_name='九张图对应数值表.xlsx')
 (out/'audit/workbook_payload.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
 mapping=[dict(figure=t['figure'],panels=t['panels'],recommendation=t['recommendation'],table=t['csv'],rows=len(t['rows'])) for t in tables]
 writecsv(out/'reports/图转表建议.csv',list(mapping[0]),[list(r.values()) for r in mapping])
 md=['# 九张图的表格化建议','',GENERAL,'','| 原图 | 建议 | 对应表格 |','|---|---|---|']
 md += [f'| {t["figure"]} {t["panels"]} | {t["recommendation"]} | {t["csv"]} |' for t in tables]
 md += ['','共1个Excel工作簿（索引＋9个数据表）和9个CSV。所有数值直接来自已采用图形的Source Data，无模型重跑。',
  '保留Fig01流程、Fig04原始分布和Fig05 PCoA，避免丢失图形所表达的结构。',
  '原Supp01是完整资料旧图，179人/10例；后续成人MI为224人/12例，分别保存。','CSV中的比例和R²为0–1数值；Excel与HTML按百分比显示。',
  '原图全部保留。NE为不可估，不用0代替；所有P为原始P。']
 (out/'reports/图转表说明_zh.md').write_text('\n'.join(md),encoding='utf-8')
 print('Prepared',len(tables),'tables;',sum(len(t['rows']) for t in tables),'source-result rows')

def verify(out):
 payload=json.loads((out/'audit/workbook_payload.json').read_text(encoding='utf-8'))
 book=out/'tables'/payload['workbook_name'];ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
 checks=[]
 def check(name,ok):
  if not ok:raise AssertionError(name)
  checks.append([name,True])
 with zipfile.ZipFile(book) as z:
  wb=ET.fromstring(z.read('xl/workbook.xml'));sheets=wb.find('m:sheets',ns)
  check('Index plus nine result worksheets',len(sheets)==10)
  strings=[]
  if 'xl/sharedStrings.xml' in z.namelist():
   strings=[''.join(i.itertext()) for i in ET.fromstring(z.read('xl/sharedStrings.xml'))]
  def val(c):
   v=c.find('m:v',ns);s=c.find('m:is',ns);typ=c.get('t')
   if typ=='inlineStr':return ''.join(s.itertext()) if s is not None else ''
   if v is None:return None
   if typ=='s':return strings[int(v.text)]
   if typ=='b':return v.text=='1'
   if typ=='str':return v.text or ''
   return float(v.text) if v.text else None
  def column(j):
   text=''
   while j:j,r=divmod(j-1,26);text=chr(65+r)+text
   return text
  for i,t in enumerate(payload['tables'],2):
   sheet=ET.fromstring(z.read(f'xl/worksheets/sheet{i}.xml'));cells={c.get('r'):c for c in sheet.findall('.//m:sheetData/m:row/m:c',ns)}
   for r,values in enumerate([t['headers']]+t['rows'],5):
    for j,expected in enumerate(values,1):
     ref=column(j)+str(r);actual=val(cells[ref]) if ref in cells else None
     if expected is None:ok=actual in (None,'')
     elif isinstance(expected,(float,int)):ok=isinstance(actual,(float,int)) and math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-14)
     else:ok=actual==expected
     check(t['sheet']+'!'+ref+' agrees with figure source',ok)
   check(t['sheet']+' no Excel error cells',not any(c.get('t')=='e' for c in cells.values()))
  check('All ten worksheet previews created',len(list((out/'audit/previews').glob('*.png')))==10)
 for r in json.loads((out/'audit/source_fingerprints.json').read_text(encoding='utf-8')):check('Source unchanged '+Path(r['path']).name,md5(Path(r['path']))==r['md5'])
 writecsv(out/'audit/validation_checks.csv',['check','passed'],checks)
 print('PASS:',len(checks),'saved-workbook/source checks')

def display(v,fmt):
 if v is None:return '—'
 if fmt=='percent1':return f'{v:.1%}'
 if fmt=='percent4':return f'{v:.4%}'
 if fmt=='p':return f'{v:.2e}' if v<.0001 else f'{v:.6f}'
 if fmt=='count':return f'{v:,}'
 if fmt in ('one','three','four'):return format(v,'.'+str({'one':1,'three':3,'four':4}[fmt])+'f')
 return str(v)

def publish(out):
 payload=json.loads((out/'audit/workbook_payload.json').read_text(encoding='utf-8'))
 target=ROOT/'Results_表格与图片/Results_重点表格/tables';target.mkdir(exist_ok=True)
 if not target.resolve().is_relative_to(ROOT.resolve()):raise ValueError('Unsafe output target')
 copy_files=[out/'tables'/payload['workbook_name']]+[out/'tables'/t['csv'] for t in payload['tables']]
 rows=[]
 for src in copy_files:
  dest=target/src.name
  if dest.exists() and md5(dest)!=md5(src):
   archive=out/'audit/previous_viewing_copies';archive.mkdir(exist_ok=True);shutil.copy2(dest,archive/dest.name)
  shutil.copy2(src,dest)
  if md5(src)!=md5(dest):raise ValueError('Copy mismatch')
  rows.append([src.name,str(src),md5(src)])
 writecsv(target/'file_manifest.csv',['file','source','md5'],rows)
 page=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>九张图对应表格</title>',
 '<style>body{font:15px/1.7 Arial,"Microsoft YaHei";margin:30px;color:#243746}table{border-collapse:collapse;white-space:nowrap}th{background:#79CAFB;color:#17202a}th,td{border-bottom:1px solid #ddd;padding:7px 12px;text-align:right}th:first-child,td:first-child{text-align:left}.scroll{overflow:auto}p{max-width:1100px}a{color:#294fa8}</style>',
 '<h1>九张Nature图对应数值表</h1><p>'+GENERAL+'</p>',f'<p><a href="{payload["workbook_name"]}">下载Excel工作簿（索引＋9个数据表）</a></p>',
 '<p>优先转表：Fig02、Fig03、Supp01–04。部分转表：Fig01 b–d、Fig04 b、Fig05 b。流程、分布和PCoA图建议保留。</p>']
 for t in payload['tables']:
  page += ['<h2>'+html.escape(t['title'])+'</h2>','<p>'+t['figure']+' '+t['panels']+' · '+t['recommendation']+'</p>',
   '<div class="scroll"><table><thead><tr>'+''.join('<th>'+html.escape(h)+'</th>' for h in t['headers'])+'</tr></thead><tbody>']
  page += ['<tr>'+''.join('<td>'+html.escape(display(v,f))+'</td>' for v,f in zip(r,t['formats']))+'</tr>' for r in t['rows']]
  page += ['</tbody></table></div>','<p>'+'<br>'.join(html.escape(n) for n in t['notes'])+'</p>',f'<p><a href="{t["csv"]}">下载CSV</a></p>']
 page += ['</html>'];(target/'index.html').write_text('\n'.join(page),encoding='utf-8')
 shutil.copy2(out/'reports/图转表说明_zh.md',target/'README_zh.md')
 main=ROOT/'Results_表格与图片/Results_重点表格/index.html';text=main.read_text(encoding='utf-8-sig')
 if 'tables/index.html' not in text:
  text=text.replace('</ul>','<li><a href="tables/index.html">九张Nature图转表：Excel与CSV</a></li></ul>')
  main.write_text(text,encoding='utf-8')
 print('Published requested folder:',target)

if __name__=='__main__':
 mode=sys.argv[1];out=Path(sys.argv[2]);{'prepare':prepare,'verify':verify,'publish':publish}[mode](out)
