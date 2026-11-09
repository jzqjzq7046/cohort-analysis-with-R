"""English journal tables from adopted aggregate results, without statistical refits."""
from pathlib import Path
import csv, hashlib, html, json, math, re, shutil, sys, zipfile
import xml.etree.ElementTree as ET
from prepare_figure_tables import current, read, number, md5, writecsv, PAIR, SOURCE_NAMES
from consolidate_key_table_folder import rebuild_index

ROOT=Path(__file__).resolve().parents[1]
VIEW=ROOT/'Results_表格与图片/Results_重点表格/tables/English'
REVIEW=ROOT/'literature_reviews/20261020_english_tables'
SCENARIOS={
 'core_base':'Technical reference', 'core_adjusted':'Primary adjustment',
 'adult_core':'Adult core', 'adult_lifestyle':'Adult lifestyle', 'adult_BMI':'Adult lifestyle + BMI',
 'bowel_base':'Bowel-data reference', 'bowel_adjusted':'Bowel symptoms + probiotics',
 'unperturbed_core':'No recent perturbation', 'DNA_base':'DNA-data reference', 'core_DNA':'Additional DNA adjustment'}
GENERAL='Baseline eligibility: 319 primary families with at least one observed member in each of G0, G1 and G2 (1,275 participants). Minors are retained in the primary analyses. All P values are unadjusted for multiple testing.'
CLINICAL=[
 'The exposure is at least one observed older-generation member with ultrasound-detected hepatic steatosis. The reference is no observed older-generation case. Each younger-generation participant contributes one row to each comparison.',
 'Primary adjusted models were adjusted for the younger participant\'s age and sex; the number, mean age and proportion of women among observed older-generation members; survey centre and urbanicity. Covariates with no variation were omitted. All models, including unadjusted models, used logistic GEE with family-clustered robust standard errors.',
 'Primary family labels identify shared family membership. Direct parent-child or grandparent-grandchild relationships and current cohabitation have not been verified.'
]
PAIR_NOTES=[
 'Primary adjustment uses both participants: age (3-d.f. natural splines), sex, survey centre, library batch, sequencing pool, log10 input-matrix total, sampling season (sine/cosine), log1p recorded collection-to-freezing interval, and recent antibiotic use. Indicators for the same centre, batch and pool are included; aliased columns are removed without changing the identifiable adjustment space.',
 'IF, within family; BF, between families. Lower Bray–Curtis dissimilarity indicates greater compositional similarity. The adjusted contrast is IF minus BF. Pair counts are not independent sample sizes because participants recur across pairs.',
 'Two-sided P values use 9,999 studentized Freedman–Lane node-residual permutations within generation-role, sex and 5-year age strata. The add-one minimum P is 0.0001. No confidence interval treating repeated pairs as independent is reported. Compositional similarity does not establish strain transmission.'
]
DIV_NOTES=[
 'Primary adjustment: age (generation-specific 3-d.f. natural splines), sex, generation and required interactions, survey centre, library batch, sequencing pool, log10 input-matrix total, sampling season (sine/cosine), log1p recorded collection-to-freezing interval, log1p absolute stool-to-ultrasound interval, and recent antibiotic use.',
 'Input-matrix total is the total of the supplied input abundance/count matrix; it has not been verified as raw sequencing read depth.'
]
SENS_NOTES=[
 'Technical reference and primary adjustment use the same primary complete-case sample. The technical reference already controls age, sex, library batch, sequencing pool and log10 input-matrix total; it is not an unadjusted model.',
 'Adult core, lifestyle and BMI models share one adult complete-case sample per analysis. Lifestyle adds alcohol, smoking, diet score, adult physical activity and log10 total energy intake; BMI is added last.',
 'Bowel-data reference and bowel-symptom models share one sample; additions are constipation, loose stools and probiotic/lactic-acid preparation use. DNA-data reference and DNA-adjusted models share one sample; log10 DNA concentration is added.',
 'No recent perturbation excludes recent antibiotic use or gastroenteritis. Sample-restricted estimates cannot be compared with full-sample estimates as a pure effect of covariate adjustment.'
]
REFERENCES=[
 dict(id='R1',citation='Long MT et al. Liver International 39, 740–747 (2019).',title='Parental non-alcoholic fatty liver disease increases risk of non-alcoholic fatty liver disease in offspring.',doi='10.1111/liv.13956',url='https://pmc.ncbi.nlm.nih.gov/articles/PMC6758911/',location='Tables 1 and 2',observed='Characteristics by exposure group; crude and successive adjusted models shown as rows with OR (95% CI) and P.',adopted='Table 1 characteristic-by-group layout; Table 2 crude and primary models in successive rows. Our covariate set is retained, not replaced by theirs.'),
 dict(id='R2',citation='Caussy C et al. Journal of Clinical Investigation 127, 2697–2704 (2017).',title='Nonalcoholic fatty liver disease with cirrhosis increases familial risk for advanced fibrosis.',doi='10.1172/JCI93465',url='https://www.jci.org/articles/view/93465',location='Tables 1 and 2; Figure 3',observed='Clinical characteristics are separated by study group; crude and adjusted associations include 95% CIs.',adopted='Keep descriptive characteristics separate from model effects; define comparison and denominator. This study does not verify our pedigree labels.'),
 dict(id='R3',citation='Gacesa R et al. Nature 604, 732–739 (2022).',title='Environmental factors shaping the gut microbiome in a Dutch population.',doi='10.1038/s41586-022-04567-7',url='https://www.nature.com/articles/s41586-022-04567-7',location='Supplementary workbook: S2A, S2B and S3A–S3C',observed='Metadata tables include available and missing counts; association tables include N, effect size, R2, F statistic, P and FDR as distinct fields.',adopted='Retain N, effect and P with explicit missingness; retain F and R2 for composition. Per user instruction we report raw P, not the paper’s FDR.'),
 dict(id='R4',citation='Valles-Colomer M et al. Nature Microbiology 7, 87–96 (2022).',title='Variation and transmission of the human gut microbiota across multiple familial generations.',doi='10.1038/s41564-021-01021-8',url='https://www.nature.com/articles/s41564-021-01021-8',location='Supplementary workbook: Legends and Table_S2',observed='The workbook separates metadata, community-variation analyses and strain analyses; Table S2 reports independent/cumulative dbRDA explanatory power and significance.',adopted='Keep Bray–Curtis compositional similarity separate from strain transmission; show explained variation with permutation P and stated sample scope. Their multiplicity correction is not adopted.'),
 dict(id='R5',citation='Lee G et al. Nature Communications 11, 4982 (2020).',title='Distinct signatures of gut microbiome and metabolites associated with significant fibrosis in non-obese NAFLD.',doi='10.1038/s41467-020-18754-5',url='https://www.nature.com/articles/s41467-020-18754-5',location='Supplementary Tables 1, 2 and 6 (PDF pp. 25, 27 and 33)',observed='Characteristics are rows; clinical subgroups and their sample sizes are columns; units and distribution summaries are defined in notes.',adopted='Table 1 uses generation columns, explicit units, median [Q1, Q3] and n/N (%). We preserve our existing summaries rather than switching to their mean/SD.'),
]

def table(code,sheet,title,figure,headers,rows,formats,widths,notes,sources,ci=None):
 assert len(headers)==len(formats)==len(widths)
 assert all(len(r)==len(headers) for r in rows)
 return dict(code=code,sheet=sheet,title=title,figure=figure,headers=headers,rows=rows,formats=formats,widths=widths,notes=notes,sources=sources,ci=ci,csv=code+'_'+sheet.replace(' ','_')+'.csv')

def prepare(out):
 src=current();sd=src/'source_data';d={k:read(sd/(k+'.csv')) for k in SOURCE_NAMES};tables=[]
 mi_runs=[p for p in (ROOT/'outputs/13_Clinical_MI/current').iterdir() if p.is_dir()]
 assert len(mi_runs)==1, 'Exactly one adopted adult MI run is required'
 mi_run=mi_runs[0]
 assert (mi_run/'run_status.txt').read_text(encoding='utf-8-sig').splitlines()[0]=='COMPLETED'
 mi_source=mi_run/'tables/pooled_adult_models.csv'
 mi=read(mi_source)
 assert len(mi)==9 and len({(r['comparison'],r['model']) for r in mi})==9
 # A compact manuscript Table 1, using existing descriptive summaries only.
 baseline=list((ROOT/'outputs/01_Table1/current').glob('*/tables/Table1_by_generation.csv'))
 assert len(baseline)==1
 b={int(r['row_id']):r for r in read(baseline[0])}
 ids=[1,7,9,4,5,15,17,18,19,20,21,24,25,26,27,28,29,30,31,34,35,37,38]
 labels={7:'Female sex, n/N (%)',9:'Age <18 years, n/N (%)',4:'BMI, kg/m²',5:'Waist circumference, cm',
 15:'Ultrasound-detected hepatic steatosis, n/N (%)',17:'Mild hepatic fat accumulation, n/N (%)',18:'Clear hepatic fat accumulation, n/N (%)',
 19:'Systolic blood pressure, mmHg',20:'Diastolic blood pressure, mmHg',34:'Alcohol intake, servings/day',35:'Diet score (0–27)',37:'Daily tobacco use, n/N (%)',38:'Antibiotics within 3 months before stool collection, n/N (%)'}
 rows=[['Participants, n',1275,468,319,488]]
 rows += [[labels.get(i,b[i]['characteristic'])]+[b[i][g] for g in ['Overall','G0','G1','G2']] for i in ids]
 tables.append(table('Table1','Cohort','Clinical characteristics and hepatic steatosis by generation','Fig01',
  ['Characteristic','Overall (N = 1,275)','G0 (N = 468)','G1 (N = 319)','G2 (N = 488)'],rows,
  ['text']*5,[370,235,235,235,235],
  ['Unadjusted descriptive summaries. Continuous variables are median [Q1, Q3]; available N is retained in each cell. Categorical variables are n/nonmissing N (%). No between-generation significance tests were added.',
   'G0, oldest; G1, middle; G2, youngest generation. Among G2 participants, steatosis was observed in 5/264 minors (1.9%) and 12/224 adults (5.4%). These groups are subsets, not additional participants.',
   'BMI, body mass index; ALT, alanine aminotransferase; AST, aspartate aminotransferase; CRP, C-reactive protein; HDL/LDL, high-/low-density lipoprotein. Steatosis denotes the recorded ultrasound outcome, not a newly adjudicated MASLD diagnosis.',
   'Formatting references: R1 Tables 1–2; R2 Tables 1–2; R5 Supplementary Tables 1–2. The original full descriptive table is retained in module 01.'],['01_Table1/Table1_by_generation.csv','proportions.csv']))
 prev={(r['comparison'],r['older']):r for r in d['clinical_prevalence_same_sample']}
 rows=[]
 for pair in PAIR:
  for model,label in [('R0_crude_same_sample','Unadjusted'),('R1_review_primary','Primary adjusted')]:
   r=next(r for r in d['clinical'] if r['comparison']==pair and r['model']==model)
   a,c=prev[(pair,'No older case')],prev[(pair,'Older case')]
   def ratio(x):return f"{x['cases']}/{x['n']} ({100*number(x,'cases')/number(x,'n'):.1f}%)"
   rows.append([pair,label,number(r,'n'),ratio(a),ratio(c)]+[number(r,k) for k in ['OR','CI_low','CI_high','p_value']])
 tables.append(table('Table2','Clinical association','Cross-generation association of hepatic steatosis','Fig02',
  ['Comparison','Model','Younger participants, N','No older case: cases/N (%)','≥1 older case: cases/N (%)','OR','95% CI lower','95% CI upper','P value'],rows,
  ['text','text','count','text','text','two','two','two','p'],[115,175,125,180,180,80,110,110,110],
  CLINICAL+['Unadjusted and primary adjusted estimates use exactly the same complete-case sample within each comparison. The outcome is steatosis in the younger generation. Unadjusted models still account for family clustering.',
  'Presentation follows R1 Table 2. All-age estimates must not be equated with adult-only MI estimates.'],['clinical.csv','clinical_prevalence_same_sample.csv'],ci=[5,6,7]))
 rows=[]
 for pair in PAIR:
  r=next(r for r in d['distance'] if r['comparison']==pair and r['scenario']=='core_adjusted')
  rows.append([pair]+[number(r,k) for k in ['n_people','families','within_pairs','between_pairs','observed_within_mean','observed_between_mean','adjusted_distance_difference','p_raw_two_sided']])
 tables.append(table('Table3','Family similarity','Within-family and between-family microbiome dissimilarity','Fig03',
  ['Comparison','Participants, N','Families, N','IF pairs, n','BF pairs, n','Mean IF distance','Mean BF distance','Adjusted IF − BF','P value'],rows,
  ['text']+['count']*4+['four']*3+['p'],[115,115,100,110,135,130,130,150,110],PAIR_NOTES+['Original IF/BF means are descriptive. Presentation separates sample counts, effect estimates and P as in R3; interpretation follows the community/strain distinction in R4.'],['distance.csv','pair.csv']))
 rows=[]
 for g in ['G0','G1','G2']:
  r=next(r for r in d['shannon'] if r['generation']==g and r['scenario']=='core_adjusted')
  rows.append([g]+[number(r,k) for k in ['n','cases','estimate','CI_low','CI_high','p_raw']])
 hp=number(next(r for r in d['heterogeneity'] if r['scenario']=='core_adjusted'),'p_raw')
 tables.append(table('Table4','Shannon','Association of hepatic steatosis with Shannon diversity','Fig04b',
  ['Generation','Participants, N','Steatosis cases, n','Adjusted difference','95% CI lower','95% CI upper','P value'],rows,
  ['text','count','count','three','three','three','p'],[115,120,135,150,120,120,110],
  DIV_NOTES+['Gaussian identity-link GEE with family-clustered robust 95% CIs. Differences are steatosis minus no steatosis. The primary sample comprises 1,267 participants and 232 cases.',f'Generation-by-steatosis heterogeneity P = {hp:.4f}. Different point estimates alone do not establish a generation difference.','Keep the original raw-distribution figure; this table reports model estimates. Presentation separates N, effect and P as in R3.'],['shannon.csv','heterogeneity.csv'],ci=[3,4,5]))
 def composition_rows(scenarios):
  rows=[]
  for scenario in scenarios:
   a=next(r for r in d['composition'] if r['scenario']==scenario and r['scheme']=='family')
   c=next(r for r in d['composition'] if r['scenario']==scenario and r['scheme']=='family_generation')
   for k in ['n','cases','families','pseudo_F','R2_total_corrected']:assert number(a,k)==number(c,k)
   rows.append([SCENARIOS[scenario]]+[number(a,k) for k in ['n','cases','families','pseudo_F','R2_total_corrected','p_raw']]+[number(c,'p_raw')])
  return rows
 comp_headers=['Adjustment','Participants, N','Steatosis cases, n','Families, N','Pseudo-F','Explained variation, %','P: within family','P: within family and generation']
 comp_formats=['text','count','count','count','three','percent4','p','p']
 comp_widths=[275,115,135,100,105,145,130,170]
 comp_notes=DIV_NOTES+['Lingoes-corrected Bray–Curtis partial dbRDA additionally conditions on family. R² is the steatosis-associated fraction of total corrected variation, displayed as a percentage (0–1 fraction in CSV).',
  'The two permutation restrictions use the same fitted model; N, pseudo-F and R² therefore appear once. P values are upper-tail permutation P from 9,999 permutations. They are not two independent fitted effects.',
  'PCoA coordinates remain unadjusted descriptive coordinates and are retained as a figure. Showing explained variation with significance follows R4 Table S2 and the distinct effect/test fields in R3.']
 tables.append(table('Table5','Composition','Conditional association with overall microbiome composition','Fig05b',comp_headers,composition_rows(['core_adjusted']),comp_formats,comp_widths,comp_notes,['composition.csv']))
 models={'R2_adult_core_same_sample':'Core','R3_adult_lifestyle':'Core + lifestyle','R4_adult_lifestyle_BMI':'Core + lifestyle + BMI'}
 rows=[]
 for pair in PAIR:
  for model,label in models.items():
   r=next(r for r in mi if r['comparison']==pair and r['model']==model)
   assert number(r,'m_total')==50
   if r['status']=='POOLED':
    assert number(r,'m_estimable')==50 and all(number(r,k) is not None for k in ['OR','CI_low','CI_high','p_raw'])
    status='Pooled (50/50)'
   else:
    assert r['status']=='NOT_POOLABLE' and r['reason']=='SEPARATION'
    assert all(number(r,k) is None for k in ['OR','CI_low','CI_high','p_raw'])
    status='Not poolable: separation'
   rows.append([pair,label]+[number(r,k) for k in ['n','cases','OR','CI_low','CI_high','p_raw']]+[status])
 tables.append(table('TableS1','Adult clinical MI','Adult clinical extensions after multiple imputation','Supp05_Adult_MI_comparison',
  ['Comparison','Model','Participants, N','Steatosis cases, n','OR','95% CI lower','95% CI upper','P value','Estimation status'],rows,
  ['text','text','count','count','two','two','two','p','text'],[115,205,115,130,85,110,110,110,225],
  CLINICAL+['Adult-only supplementary analyses. Core uses the primary clinical adjustment set. Lifestyle adds alcohol, smoking, diet score, adult physical activity and log10 total energy intake; the last model additionally adjusts BMI. All three models retain the same eligible adult sample within each comparison.',
  'Missing adult covariates were multiply imputed using family-multilevel MICE (50 completed datasets). Logistic GEE was fitted in each dataset; log ORs and family-robust variances were combined using Rubin rules. Participants and cases are unique people, not totals across imputations. Core rows use the same adult sample as the extension models.',
  'Only models estimable in all 50 datasets are pooled. G1-G2 with BMI was estimable in 21/50 datasets and G0-G2 with BMI in 0/50; separation prevents valid pooling. Blank numeric cells denote unavailable estimates, not zero.',
  'This version replaces the original complete-case Table S1. Original complete-case estimates remain in archived module-15 outputs and module 13 complete_case_vs_MI.csv. The all-age primary results in Table 2 are unchanged.',
  'Successive model rows follow R1 Table 2; failure status is an additional reporting requirement for these data.'],[str(mi_source.relative_to(ROOT))],ci=[4,5,6]))
 rows=[]
 for scenario,label in SCENARIOS.items():
  for pair in PAIR:
   r=next(r for r in d['distance'] if r['comparison']==pair and r['scenario']==scenario)
   rows.append([label,pair]+[number(r,k) for k in ['n_people','families','adjusted_distance_difference','p_raw_two_sided']])
 tables.append(table('TableS2','Similarity sensitivity','Sensitivity analyses of family microbiome similarity','Supp02',
  ['Analysis','Comparison','Participants, N','Families, N','Adjusted IF − BF','P value'],rows,
  ['text','text','count','count','four','p'],[275,115,120,110,160,120],SENS_NOTES+PAIR_NOTES,['distance.csv']))
 rows=[]
 for scenario,label in SCENARIOS.items():
  for g in ['G0','G1','G2']:
   r=next(r for r in d['shannon'] if r['generation']==g and r['scenario']==scenario)
   rows.append([label,g]+[number(r,k) for k in ['n','cases','estimate','CI_low','CI_high','p_raw']])
 tables.append(table('TableS3','Shannon sensitivity','Sensitivity analyses of Shannon diversity associations','Supp03',
  ['Analysis','Generation','Participants, N','Steatosis cases, n','Adjusted difference','95% CI lower','95% CI upper','P value'],rows,
  ['text','text','count','count','three','three','three','p'],[275,105,115,130,145,120,120,110],
  SENS_NOTES+DIV_NOTES+['Differences are steatosis minus no steatosis, with robust family-clustered 95% CIs and two-sided raw P from Gaussian GEE. Adult rows are the original complete-case analyses, not MI.'],['shannon.csv'],ci=[4,5,6]))
 tables.append(table('TableS4','Composition sensitivity','Sensitivity analyses of overall microbiome composition','Supp04',
  comp_headers,composition_rows(SCENARIOS),comp_formats,comp_widths,SENS_NOTES+comp_notes,['composition.csv']))
 assert len(tables)==9
 for t in tables:
  assert not re.search('[\u4e00-\u9fff]',json.dumps([t['title'],t['headers'],t['rows'],t['notes']],ensure_ascii=False))
  writecsv(out/'tables'/t['csv'],t['headers'],t['rows'])
 protected=[p for p in (src/'figures').iterdir() if p.is_file()]+[sd/(k+'.csv') for k in SOURCE_NAMES]+baseline+[mi_source]
 hashes=[dict(path=str(p),md5=md5(p)) for p in protected]
 (out/'audit/source_fingerprints.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
 payload=dict(workbook_name='Hepatic_steatosis_journal_tables_EN.xlsx',source_run=src.name,general=GENERAL,tables=tables,references=REFERENCES)
 (out/'audit/workbook_payload.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
 write_literature_report(out)
 make_html(out,payload)
 print('Prepared',len(tables),'English tables;',sum(len(t['rows']) for t in tables),'rows; no model refits')

def ptext(v):
 if v is None:return '—'
 return f'{v:.2e}' if v<.001 else f'{v:.3g}'

def display(v,f):
 if v is None:return '—'
 if f=='p':return ptext(v)
 if f=='percent4':return f'{v:.4%}'
 if f=='count':return f'{v:,}'
 if f in ('two','three','four'):return f"{v:.{dict(two=2,three=3,four=4)[f]}f}"
 return str(v)

def make_html(out,payload):
 css='body{font:15px/1.6 Arial,sans-serif;color:#202020;max-width:1380px;margin:36px auto;padding:0 25px}h1{font-size:24px}h2{font-size:18px;margin-top:38px}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;border-top:2px solid #222;border-bottom:2px solid #222;font-size:13px}th{border-bottom:1px solid #222;font-weight:600;padding:9px 12px;text-align:center}td{padding:7px 12px;vertical-align:top;text-align:right;white-space:nowrap}td.text{text-align:left;white-space:normal}tr.group td{border-top:1px solid #aaa}.notes{font-size:12px;max-width:1250px;margin-top:12px}a{color:#245882}@media print{body{margin:0;padding:0;font-size:10pt}.scroll{overflow:visible}article{break-before:page}article:first-of-type{break-before:auto}nav,.downloads{display:none}thead{display:table-header-group}tr{break-inside:avoid}table{font-size:9pt}td,th{padding:5px}}'
 head='<!doctype html><html lang="en"><meta charset="utf-8"><title>Hepatic steatosis: English tables</title><style>'+css+'</style><body>'
 page=[head,'<h1>Hepatic steatosis across three generations</h1><p>'+GENERAL+'</p>',
  '<p>Five main tables and four supplementary tables. These are presentation revisions of adopted results; no statistical model was refitted.</p>',
  '<p class="downloads"><a href="'+payload['workbook_name']+'">Download editable Excel workbook</a> · <a href="Table_presentation_literature_EN.md">Literature and presentation decisions</a></p>',
  '<nav>'+ ' · '.join(f'<a href="#{t["code"]}">{t["code"]}</a>' for t in payload['tables'])+'</nav>']
 for t in payload['tables']:
  headers=t['headers'][:];formats=t['formats'][:];rr=[r[:] for r in t['rows']]
  if t['ci']:
   effect,low,high=t['ci'];headers[effect]=('OR' if headers[effect]=='OR' else 'Adjusted difference')+' (95% CI)'
   for r in rr:
    r[effect]=('NE' if r[effect] is None else display(r[effect],formats[effect])+' ('+display(r[low],formats[low])+', '+display(r[high],formats[high])+')')
    del r[high];del r[low]
   formats[effect]='text'
   for col in [high,low]:del headers[col];del formats[col]
  part=['<article id="'+t['code']+'"><h2>'+html.escape(t['code'].replace('TableS','Supplementary Table ').replace('Table','Table '))+'. '+html.escape(t['title'])+'</h2>',
   '<div class="scroll"><table><thead><tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in headers)+'</tr></thead><tbody>']
  last=None
  for r in rr:
   cls=' class="group"' if last is not None and r[0]!=last else '';last=r[0]
   part.append('<tr'+cls+'>'+''.join('<td'+(' class="text"' if f=='text' else '')+'>'+html.escape(display(v,f))+'</td>' for v,f in zip(r,formats))+'</tr>')
  part+=['</tbody></table></div><div class="notes">'+''.join('<p>'+html.escape(n)+'</p>' for n in t['notes'])+'</div>',
   '<p class="notes">Figure correspondence: '+t['figure']+'. Sources: '+', '.join(t['sources'])+'.</p>',
   '<p class="downloads"><a href="'+t['csv']+'">Numeric CSV</a></p></article>']
  (out/'tables'/(t['code']+'.html')).write_text(head+'\n'.join(part)+'</body></html>',encoding='utf-8')
  page+=part
 page+=['<h2>Presentation references</h2><ol>'+''.join('<li>'+html.escape(r['citation'])+' <a href="'+r['url']+'">'+html.escape(r['location'])+'</a></li>' for r in REFERENCES)+'</ol></body></html>']
 (out/'tables/index.html').write_text('\n'.join(page),encoding='utf-8')

def write_literature_report(out):
 text=['# Literature-informed table presentation','',
 'Targeted review performed on 20 October 2026. This is a presentation review, not a systematic review or a new statistical analysis. Publisher/PMC articles and publisher supplements were used.','',
 '## Evidence and adaptation','',
 '| Reference | Location inspected | Reporting feature | Applied here |','|---|---|---|---|']
 for r in REFERENCES:text.append(f"| [{r['id']}: {r['citation']}]({r['url']}) | {r['location']} | {r['observed']} | {r['adopted']} |")
 text += ['', '## Presentation decisions','',
  '- Tables 1–5 correspond to the five main figures; Tables S1–S4 correspond to the four supplementary figures. Table 1 also uses existing clinical summaries from module 01. The full original Table 1 is retained.',
  '- HTML tables use a top rule, a header rule and a bottom rule, without vertical borders. This is our restrained manuscript presentation choice, not a claim that every cited journal mandates one universal format.',
  '- HTML combines each point estimate and 95% CI in one field. The editable Excel and numeric CSV retain point/lower/upper values as separate numeric columns so precision is not lost.',
  '- Baseline variables are rows, generations are columns. Existing medians and quartiles are preserved; n/nonmissing N (%) makes differing denominators explicit. No new baseline P values are computed.',
  '- Clinical models are successive rows within each generation comparison, as in Long Table 2. Crude and primary estimates use the same sample.',
  '- Overall-composition tables show pseudo-F, R² as percent and both restricted-permutation P values alongside a single set of N/case/family counts. Identical fits are not duplicated.',
  '- All P values remain raw per the user’s instruction. Papers using FDR or adjusted P are design references only; their multiplicity procedures are not silently imported.',
  '- Table S1 now reports adult multiple-imputation extensions from module 13. Original complete-case results remain archived and in complete_case_vs_MI.csv. Non-poolable models retain the separation reason; missing estimates are never zero.',
  '- All tables and model notes are English. Original data, images and statistical results remain unchanged.','',
  '## Retrieval record','',
  'Long and Caussy article table locations were checked through PMC/publisher pages. Gacesa and Valles-Colomer supplementary XLSX files and Lee supplementary PDF were downloaded from the publisher and inspected locally. Supplement headers and PDF table pages are retained under literature_reviews/20261020_english_tables. Some direct web-tool supplement requests failed on file type/size; local publisher retrieval succeeded.','',
  '## References','']
 for r in REFERENCES:text += [r['id']+'. '+r['citation']+' '+r['title']+' https://doi.org/'+r['doi'],'']
 report='\n'.join(text)
 (out/'reports/Table_presentation_literature_EN.md').write_text(report,encoding='utf-8')
 (REVIEW/'Table_presentation_literature_EN.md').write_text(report,encoding='utf-8')
 (out/'tables/Table_presentation_literature_EN.md').write_text(report,encoding='utf-8')
 (out/'reports/README_EN.md').write_text('# English journal tables\n\n'+GENERAL+'\n\nOpen tables/index.html for manuscript-style tables; Excel/CSV preserve numeric values. Independent entry: scripts/15_export_english_journal_tables.R. Original images and raw data are preserved. Table S1 reports adult extensions after 50-fold multiple imputation. Original complete-case estimates remain in module-15 history and module-13 complete_case_vs_MI.csv. The all-age primary Table 2 is unchanged.\n\nCanonical viewing folder: Results_表格与图片/Results_重点表格/tables/English.\n',encoding='utf-8')

def verify(out):
 payload=json.loads((out/'audit/workbook_payload.json').read_text(encoding='utf-8'));checks=[]
 def check(name,ok):
  if not ok:raise AssertionError(name)
  checks.append([name,True])
 ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
 with zipfile.ZipFile(out/'tables'/payload['workbook_name']) as z:
  wb=ET.fromstring(z.read('xl/workbook.xml'));check('Index and nine English result sheets',len(wb.find('m:sheets',ns))==10)
  strings=[]
  if 'xl/sharedStrings.xml' in z.namelist():strings=[''.join(i.itertext()) for i in ET.fromstring(z.read('xl/sharedStrings.xml'))]
  def value(c):
   typ=c.get('t');v=c.find('m:v',ns)
   if typ=='inlineStr':return ''.join(c.find('m:is',ns).itertext())
   if v is None:return None
   if typ=='s':return strings[int(v.text)]
   if typ=='str':return v.text
   return float(v.text) if v.text else None
  for i,t in enumerate(payload['tables'],2):
   root=ET.fromstring(z.read(f'xl/worksheets/sheet{i}.xml'));cells={c.get('r'):c for c in root.findall('.//m:sheetData/m:row/m:c',ns)}
   for row,values in enumerate([t['headers']]+t['rows'],5):
    for j,expected in enumerate(values):
     ref=chr(65+j)+str(row);actual=value(cells[ref]) if ref in cells else None
     if expected is None:ok=actual in (None,'')
     elif isinstance(expected,(float,int)):ok=isinstance(actual,(float,int)) and math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-14)
     else:ok=actual==expected
     check(t['code']+'!'+ref+' equals source-backed payload',ok)
   check(t['code']+' has no Excel error cells',not any(c.get('t')=='e' for c in cells.values()))
 for r in json.loads((out/'audit/source_fingerprints.json').read_text(encoding='utf-8')):check('Unchanged source '+Path(r['path']).name,md5(Path(r['path']))==r['md5'])
 check('All ten worksheets rendered',len(list((out/'audit/previews').glob('*.png')))==10)
 check('Nine standalone manuscript tables',len(list((out/'tables').glob('Table*.html')))==9)
 writecsv(out/'audit/validation_checks.csv',['check','passed'],checks)
 print('PASS:',len(checks),'saved-workbook and source checks')

def publish(out):
 VIEW.mkdir(parents=True,exist_ok=True)
 assert VIEW.resolve().is_relative_to(ROOT.resolve())
 manifest=[]
 for src in (out/'tables').iterdir():
  if not src.is_file() or src.suffix.lower() not in ('.html','.xlsx','.csv','.md'):continue
  dest=VIEW/src.name
  if dest.exists() and md5(src)!=md5(dest):
   backup=out/'audit/previous_viewing_copies';backup.mkdir(exist_ok=True);shutil.copy2(dest,backup/dest.name)
  shutil.copy2(src,dest);assert md5(src)==md5(dest)
  manifest.append([src.name,str(src),md5(src)])
 writecsv(VIEW/'file_manifest.csv',['file','source','md5'],manifest)
 shutil.copy2(out/'reports/README_EN.md',VIEW/'README_EN.md')
 # Retire the superseded CSV from the managed viewing folder after the new
 # workbook and Table S1 have been copied successfully; keep its exact bytes.
 stale=VIEW/'TableS1_Adult_complete_cases.csv'
 if stale.exists():
  backup=out/'audit/previous_viewing_copies';backup.mkdir(exist_ok=True)
  shutil.copy2(stale,backup/stale.name)
  assert md5(stale)==md5(backup/stale.name)
  stale.unlink()
 rebuild_index()
 print('Published English tables:',VIEW)

if __name__=='__main__':
 command,out=sys.argv[1],Path(sys.argv[2]);globals()[command](out)
