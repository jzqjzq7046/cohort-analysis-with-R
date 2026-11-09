// Workbook authoring uses only the configured bundled @oai/artifact-tool runtime.
import fs from 'node:fs/promises';
import path from 'node:path';
import {Workbook, SpreadsheetFile} from '@oai/artifact-tool';

const [payloadPath,out]=process.argv.slice(2);
if(!payloadPath||!out)throw new Error('Usage: builder payload.json output-run');
const payload=JSON.parse(await fs.readFile(payloadPath,'utf8'));
const wb=Workbook.create();
const font='Microsoft YaHei';
const formats={text:'@',count:'#,##0',one:'0.0',three:'0.000',four:'0.0000',percent1:'0.0%',percent4:'0.0000%',p:'[<0.0001]0.00E+00;0.000000'};
const previews=path.join(out,'audit','previews');await fs.mkdir(previews,{recursive:true});
const index=wb.worksheets.add('00_图转表索引');
const review=[];
function column(n){let s='';while(n){const r=(n-1)%26;s=String.fromCharCode(65+r)+s;n=Math.floor((n-1)/26);}return s;}
function styleBase(s,last,nrows){
 s.showGridLines=false;
 s.getRange(`A1:${last}${nrows}`).format.font={name:font,size:10,color:'#17202A'};
 s.getRange(`A1:${last}${nrows}`).format.rowHeight=24;
 s.getRange(`A1:${last}${nrows}`).format.verticalAlignment='center';
 s.getRange(`A2:${last}2`).format.rowHeight=32;
 s.getRange('A2').format.font={name:font,size:14,bold:true,color:'#17202A'};
 s.getRange(`A4:${last}4`).format.borders={bottom:{style:'thin',color:'#5292F7'}};
}
function splitNote(text,width=78){
 const letters=Array.from(text),parts=[];for(let k=0;k<letters.length;k+=width)parts.push(letters.slice(k,k+width).join(''));return parts;
}
styleBase(index,'D',23);
index.getRange('A2').values=[['九张Nature图：表格化结果索引']];
index.getRange('A3').values=[['原图保留；本表忠实对应原9张图的数值版本。']];
index.getRange('A5:D14').values=[['原图','转表面板','展示建议','工作表'],...payload.tables.map(t=>[t.figure,t.panels,t.recommendation,t.sheet])];
const iw=[320,115,430,280];iw.forEach((w,i)=>index.getRange(`${column(i+1)}1:${column(i+1)}23`).format.columnWidthPx=w);
index.getRange('A5:D5').format.fill='#79CAFB';index.getRange('A5:D5').format.font.bold=true;
index.getRange('A5:D5').format.horizontalAlignment='center';
index.getRange('A17').values=[['优先整体转表：Fig02、Fig03、Supp01–04。流程、原始分布与PCoA建议继续保留。']];
index.getRange('A18').values=[['版本：原成人补充图为完整资料样本（涉及G2：179人/10例），不混入之后的成人MI。']];
index.getRange('A19').values=[['所有P为原始P。空数值与不可估状态并列展示；不会用0代替NE。']];
index.getRange('A20').values=[['比例/R²：Excel与HTML显示百分比，CSV保留0–1比例。每个结果表对应一个原图。']];
review.push({sheet:index.name,range:'A1:D21',preview:'00_图转表索引.png'});

for(const t of payload.tables){
 const s=wb.worksheets.add(t.sheet);const cols=t.headers.length,last=column(cols),end=5+t.rows.length;
 const notes=t.notes.flatMap(n=>splitNote(n,Math.max(55,Math.floor(t.widths.reduce((a,b)=>a+b,0)/11))));
 const finalRow=end+notes.length+6;styleBase(s,last,finalRow);
 s.getRange('A2').values=[[t.title]];
 s.getRange('A3').values=[[`${t.figure} ${t.panels}；${t.recommendation}`]];
 s.getRange(`A5:${last}${end}`).values=[t.headers,...t.rows];
 s.getRange(`A5:${last}5`).format={fill:'#79CAFB',font:{name:font,size:10,bold:true,color:'#17202A'},horizontalAlignment:'center',verticalAlignment:'center',wrapText:true,rowHeight:39};
 s.getRange(`A6:${last}${end}`).format.horizontalAlignment='right';
 t.widths.forEach((w,j)=>{
  const c=column(j+1);s.getRange(`${c}1:${c}${finalRow}`).format.columnWidthPx=w;
  s.getRange(`${c}6:${c}${end}`).setNumberFormat(formats[t.formats[j]]);
  if(t.formats[j]==='text')s.getRange(`${c}6:${c}${end}`).format.horizontalAlignment='left';
 });
 for(let r=1;r<t.rows.length;r++)if(t.rows[r][0]!==t.rows[r-1][0]){
  s.getRange(`A${r+6}:${last}${r+6}`).format.borders={top:{style:'thin',color:'#D9E2EC'}};
 }
 const status=t.headers.indexOf('模型状态');
 if(status>=0){
  const c=column(status+1);
  s.getRange(`${c}6:${c}${end}`).conditionalFormats.add('containsText',{text:'不可估',format:{fill:'#FBEB66'}});
 }
 let n=end+2;
 s.getRange(`A${n}`).values=[['说明']];s.getRange(`A${n}`).format.font.bold=true;
 for(const line of notes)s.getRange(`A${++n}`).values=[[line]];
 s.getRange(`A${n+2}`).values=[['来源：'+t.figure+' / '+t.sources.join('；')+'（'+payload.source_run+'）']];
 if(t.rows.length>15){s.freezePanes.freezeRows(5);s.freezePanes.freezeColumns(2);}
 review.push({sheet:s.name,range:`A1:${last}${n+2}`,preview:t.sheet+'.png'});
}
wb.recalculate();
const overview=await wb.inspect({kind:'table',range:'02_脂肪肝跨代关联!A5:K8',include:'values,formulas',tableMaxRows:4,tableMaxCols:11,maxChars:3600});
await fs.writeFile(path.join(out,'audit','workbook_inspection.ndjson'),overview.ndjson,'utf8');
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:30},summary:'Formula error scan'});
await fs.writeFile(path.join(out,'audit','formula_error_scan.ndjson'),errors.ndjson,'utf8');
for(const r of review){
 const blob=await wb.render({sheetName:r.sheet,range:r.range,scale:1,format:'png'});
 await fs.writeFile(path.join(previews,r.preview),new Uint8Array(await blob.arrayBuffer()));
 console.log('Rendered:',r.sheet);
}
const output=await SpreadsheetFile.exportXlsx(wb);await output.save(path.join(out,'tables',payload.workbook_name));
console.log('Exported workbook:',payload.workbook_name);
