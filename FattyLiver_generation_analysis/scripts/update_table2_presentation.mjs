// Targeted existing-workbook text edit; numerical result cells are unchanged.
import fs from 'node:fs/promises';
import path from 'node:path';
import {FileBlob,SpreadsheetFile} from '@oai/artifact-tool';
const [phase,input,out,payloadFile]=process.argv.slice(2);
const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(input));
const s=wb.worksheets.getItem('Clinical association');
await fs.mkdir(path.join(out,'audit','previews'),{recursive:true});
if(phase==='help'){
  console.log(wb.help('workbook.render',{include:'index,examples,notes',maxChars:6000}).ndjson);
}else if(phase==='inspect'){
  const x=await wb.inspect({kind:'table',range:'Clinical association!A5:I34',include:'values,formulas',tableMaxRows:30,tableMaxCols:9,maxChars:12000});
  await fs.writeFile(path.join(out,'audit','before_table2_inspection.ndjson'),x.ndjson);
  const image=await wb.render({sheetName:s.name,range:'A1:I32',format:'png',scale:1.3,headers:false});
  await fs.writeFile(path.join(out,'audit','previews','Before_Table2.png'),new Uint8Array(await image.arrayBuffer()));
  console.log('BEFORE_RENDERED');
}else if(phase==='apply'){
  const payload=JSON.parse(await fs.readFile(payloadFile,'utf8'));
  const t=payload.tables.find(t=>t.code==='Table2');
  for(const row of [6,8,10]){
    if(s.getRange(`B${row}`).values[0][0]!=='Crude')throw Error('Unexpected original model label');
    s.getRange(`B${row}`).values=[['Unadjusted']];
  }
  const max=Math.floor(t.widths.reduce((a,b)=>a+b,0)/7.1);
  function lines(text){const result=[];let line='';for(const word of text.split(' ')){
    if(line.length+word.length+1>max){result.push(line);line=word;}else line+=(line?' ':'')+word;
  }if(line)result.push(line);return result;}
  s.getRange('A14:I40').values=Array.from({length:27},()=>Array(9).fill(null));
  let row=14;
  for(const text of [...t.notes,`Source: ${payload.source_run}; ${t.sources.join(', ')}.`]){
    for(const line of lines(text))s.getRange(`A${row++}`).values=[[line]];
  }
  wb.recalculate();
  const check=await wb.inspect({kind:'table',range:`Clinical association!A5:I${row}`,include:'values,formulas',tableMaxRows:40,tableMaxCols:9,maxChars:15000});
  await fs.writeFile(path.join(out,'audit','table2_edit_inspection.ndjson'),check.ndjson);
  const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20}});
  await fs.writeFile(path.join(out,'audit','table2_edit_error_scan.ndjson'),errors.ndjson);
  const image=await wb.render({sheetName:s.name,range:`A1:I${row}`,format:'png',scale:1.5,headers:false});
  const bytes=new Uint8Array(await image.arrayBuffer());
  await fs.mkdir(path.join(out,'figures'),{recursive:true});
  await fs.writeFile(path.join(out,'figures','Table2_with_model_notes.png'),bytes);
  await fs.writeFile(path.join(out,'audit','previews','Table2.png'),bytes);
  await (await SpreadsheetFile.exportXlsx(wb)).save(path.join(out,'tables',payload.workbook_name));
  console.log('TABLE2_EDIT_EXPORTED');
}else throw Error('Use inspect or apply');
