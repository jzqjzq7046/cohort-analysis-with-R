// Journal-table workbook. Statistical numbers are copied from adopted sources, not refitted.
import fs from 'node:fs/promises';
import path from 'node:path';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';
const [payloadPath,out]=process.argv.slice(2);
const payload=JSON.parse(await fs.readFile(payloadPath,'utf8'));
const wb=Workbook.create(),font='Arial',review=[];
const formats={text:'@',count:'#,##0',two:'0.00',three:'0.000',four:'0.0000',percent4:'0.0000%',p:'[<0.01]0.00E+00;0.0000'};
function col(n){let s='';while(n){const r=(n-1)%26;s=String.fromCharCode(65+r)+s;n=Math.floor((n-1)/26);}return s;}
function lines(text,max){
 const words=text.split(' '),result=[];let row='';
 for(const w of words){if(row.length+w.length+1>max){result.push(row);row=w;}else row+=(row?' ':'')+w;}
 if(row)result.push(row);return result;
}
function base(s,last,end){
 s.showGridLines=false;s.getRange(`A1:${last}${end}`).format={font:{name:font,size:10,color:'#202020'},rowHeight:23,verticalAlignment:'center'};
 s.getRange('A2').format.font={name:font,size:14,bold:true,color:'#202020'};
 s.getRange(`A2:${last}2`).format.rowHeight=31;
}
const index=wb.worksheets.add('Table index');base(index,'D',36);
index.getRange('A2').values=[['Hepatic steatosis across three generations']];
index.getRange('A3').values=[['English journal tables; 319 primary families. Original figures and source data are preserved.']];
index.getRange('A5:D14').values=[['Table','Content','Original figure','Scope'],...payload.tables.map(t=>[t.code,t.title,t.figure,t.code==='TableS1'?'Adult multiple imputation':t.code.startsWith('TableS')?'Sensitivity analysis':'Primary / descriptive'])];
[115,580,140,255].forEach((w,i)=>index.getRange(`${col(i+1)}1:${col(i+1)}36`).format.columnWidthPx=w);
index.getRange('A5:D5').format={font:{name:font,size:10,bold:true},horizontalAlignment:'center',rowHeight:27,borders:{top:{style:'thin',color:'#202020'},bottom:{style:'thin',color:'#202020'}}};
index.getRange('A14:D14').format.borders={bottom:{style:'thin',color:'#202020'}};
const intro=[payload.general,'HTML combines estimates with 95% CIs for manuscript display. Excel and CSV retain the point/lower/upper values as numeric columns.','Table S1 reports adult extensions after 50-fold multiple imputation. Original complete-case results are retained separately. Table 2 retains the all-age primary analysis.','All P values are raw. Technical reference models are already adjusted. Unavailable estimates are blank and accompanied by an explicit status.','Presentation references (layout only; see the accompanying literature report for exact adaptations):'];
let r=17;for(const n of intro)for(const line of lines(n,147))index.getRange(`A${r++}`).values=[[line]];
for(const ref of payload.references)for(const line of lines(ref.id+': '+ref.citation+' '+ref.location+'; DOI '+ref.doi,147))index.getRange(`A${r++}`).values=[[line]];
review.push({sheet:index.name,range:`A1:D${r}`,file:'00_Table_index.png'});
for(const t of payload.tables){
 const s=wb.worksheets.add(t.sheet),last=col(t.headers.length),end=5+t.rows.length;
 const max=Math.floor(t.widths.reduce((a,b)=>a+b,0)/7.1);
 const notes=t.notes.flatMap(n=>lines(n,max));base(s,last,end+notes.length+8);
 s.getRange('A2').values=[[t.code.replace('TableS','Supplementary Table ').replace(/^Table(\d)/,'Table $1')+'. '+t.title]];
 s.getRange('A3').values=[[t.code==='TableS1'?'Adult supplementary analysis; 50 imputations. All-age primary results remain in Table 2.':'Correspondence: '+t.figure+'. All P values are raw.']];
 s.getRange(`A5:${last}${end}`).values=[t.headers,...t.rows];
 s.getRange(`A5:${last}5`).format={font:{name:font,size:10,bold:true},horizontalAlignment:'center',verticalAlignment:'center',wrapText:true,rowHeight:43,borders:{top:{style:'thin',color:'#202020'},bottom:{style:'thin',color:'#202020'}}};
 t.widths.forEach((w,i)=>{
  const c=col(i+1);s.getRange(`${c}1:${c}${end+notes.length+8}`).format.columnWidthPx=w;
  s.getRange(`${c}6:${c}${end}`).setNumberFormat(formats[t.formats[i]]);
  s.getRange(`${c}6:${c}${end}`).format.horizontalAlignment=t.formats[i]==='text'?'left':'right';
 });
 if(t.code==='Table2')s.getRange(`D6:E${end}`).format.horizontalAlignment='right';
 if(t.code==='Table1'){
  s.getRange(`B6:E6`).setNumberFormat('#,##0');s.getRange(`B6:E${end}`).format.horizontalAlignment='right';
  s.getRange(`A6:A${end}`).format.wrapText=true;
  s.getRange(`A6:${last}${end}`).format.rowHeight=27;
  s.getRange('A29:E29').format.rowHeight=38;
 }
 for(let i=1;i<t.rows.length;i++)if(t.code!=='Table1' && t.rows[i][0]!==t.rows[i-1][0])s.getRange(`A${i+6}:${last}${i+6}`).format.borders={top:{style:'thin',color:'#C6C6C6'}};
 s.getRange(`A${end}:${last}${end}`).format.borders={bottom:{style:'thin',color:'#202020'}};
 let n=end+2;s.getRange(`A${n}`).values=[['Notes']];s.getRange(`A${n}`).format.font.bold=true;
 for(const line of notes)s.getRange(`A${++n}`).values=[[line]];
 for(const line of lines('Source: '+payload.source_run+'; '+t.sources.join(', ')+'.',max))s.getRange(`A${++n}`).values=[[line]];
 if(t.rows.length>15){s.freezePanes.freezeRows(5);s.freezePanes.freezeColumns(1);}
 review.push({sheet:s.name,range:`A1:${last}${n+1}`,file:t.code+'.png'});
}
wb.recalculate();
const inspect=await wb.inspect({kind:'table',range:'Clinical association!A5:I11',include:'values,formulas',tableMaxRows:7,tableMaxCols:9,maxChars:3600});
await fs.writeFile(path.join(out,'audit','workbook_inspection.ndjson'),inspect.ndjson);
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},summary:'Final error scan'});
await fs.writeFile(path.join(out,'audit','error_scan.ndjson'),errors.ndjson);
await fs.mkdir(path.join(out,'audit','previews'),{recursive:true});
for(const p of review){
 const blob=await wb.render({sheetName:p.sheet,range:p.range,scale:1,format:'png'});
 await fs.writeFile(path.join(out,'audit','previews',p.file),new Uint8Array(await blob.arrayBuffer()));console.log('Rendered',p.sheet);
}
const output=await SpreadsheetFile.exportXlsx(wb);await output.save(path.join(out,'tables',payload.workbook_name));
console.log('Exported English workbook');
