"""Retrieve publisher supplementary tables for a targeted table-presentation review."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import urllib.request, json, zipfile, xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT/'literature_reviews/20261020_english_tables'
SOURCES = {
 'Gacesa2022.xlsx': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-022-04567-7/MediaObjects/41586_2022_4567_MOESM4_ESM.xlsx',
 'VallesColomer2022.xlsx': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-021-01021-8/MediaObjects/41564_2021_1021_MOESM2_ESM.xlsx',
 'Lee2020_supplement.pdf': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-020-18754-5/MediaObjects/41467_2020_18754_MOESM1_ESM.pdf',
}

def get(item):
 name,url=item; dest=DEST/'sources'/name
 try:
  if not dest.exists():
   req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'})
   with urllib.request.urlopen(req,timeout=45) as response: data=response.read()
   if not data.startswith(b'PK' if name.endswith('xlsx') else b'%PDF'):raise ValueError('Unexpected response format')
   dest.write_bytes(data)
  return dict(file=name,url=url,status='retrieved',bytes=dest.stat().st_size)
 except Exception as e:return dict(file=name,url=url,status='failed',error=str(e))

def inspect_xlsx(p):
 ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}; result=[]
 with zipfile.ZipFile(p) as z:
  strings=[]
  if 'xl/sharedStrings.xml' in z.namelist():
   strings=[''.join(t.text or '' for t in r.findall('.//m:t',ns)) for r in ET.fromstring(z.read('xl/sharedStrings.xml'))]
  wb=ET.fromstring(z.read('xl/workbook.xml'))
  rels={r.get('Id'):r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
  for s in wb.find('m:sheets',ns):
   rid=s.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id');target=rels[rid]
   name=target.lstrip('/') if target.startswith('/') else 'xl/'+target
   root=ET.fromstring(z.read(name)); rows=[]
   for row in root.findall('.//m:sheetData/m:row',ns)[:9]:
    vals=[]
    for c in list(row)[:16]:
     v=c.find('m:v',ns);txt=v.text if v is not None else ''
     if c.get('t')=='s':txt=strings[int(txt)]
     if c.get('t')=='inlineStr':txt=''.join(c.itertext())
     vals.append(txt)
    rows.append(vals)
   result.append(dict(sheet=s.get('name'),first_rows=rows))
 return result

if __name__=='__main__':
 (DEST/'sources').mkdir(parents=True,exist_ok=True)
 with ThreadPoolExecutor(max_workers=3) as pool:log=list(pool.map(get,SOURCES.items()))
 (DEST/'retrieval_log.json').write_text(json.dumps(log,indent=2),encoding='utf-8')
 for r in log:
  print(r['file'],r['status'],r.get('bytes',r.get('error')))
  if r['status']=='retrieved' and r['file'].endswith('xlsx'):
   evidence=inspect_xlsx(DEST/'sources'/r['file'])
   (DEST/(r['file']+'.headers.json')).write_text(json.dumps(evidence,indent=2),encoding='utf-8')
   print('Worksheets:',len(evidence))
  elif r['status']=='retrieved':
   from pypdf import PdfReader
   reader=PdfReader(DEST/'sources'/r['file'])
   pages=[(i+1,p.extract_text()) for i,p in enumerate(reader.pages)]
   relevant=[(i,t) for i,t in pages if 'Supplementary Table' in t or 'Supplementary table' in t]
   (DEST/'Lee2020_table_text.txt').write_text('\n'.join(f'PAGE {i}\n{t}' for i,t in relevant),encoding='utf-8')
   print('Table-containing pages:',[i for i,_ in relevant])
