"""Update active entry points; originals are preserved by bootstrap backup."""
from pathlib import Path
import json
root=Path(__file__).resolve().parents[1]
backup=root/'code_history/before_two_generation_20261024'
mapping={'00':'prepare','01':'descriptive','02':'visualization','03':'descriptive','04':'descriptive','05':'prepare','06':'clinical','07':'similarity','08':'diversity','09':'prepare','10':'visualization','11':'visualization','12':'clinical_mi','13':'clinical_mi','14':'tables','15':'tables','16':'visualization','17':'taxa','18':'taxa','19':'taxa','20':'visualization','21':'taxa','22':'taxa','23':'taxa'}
manifest=[]
for old in sorted((backup/'scripts').glob('[0-9][0-9]_*.R')):
    stage=mapping.get(old.name[:2])
    if not stage:continue
    note=''
    if old.name[:2] in ['12','13']:
        note="message('Clinical C2 MI is rebuilt for the current two-generation scope. Only missing adult alcohol/BMI are imputed; old restricted-cohort MI data are not reused.')\n"
    text=f'''# Current C1/C2/B1/B2/P1/P2 entry, 2026-10-24.
# Original source retained at code_history/before_two_generation_20261024/scripts/{old.name}
# Historical method-specific filenames now delegate to the adopted protocol.
.current_root <- Sys.getenv("FATTYLIVER_ROOT", unset="E:/YFS/workspace/FattyLiver_generation_analysis")
if(.Platform$OS.type=="windows")invisible(Sys.setlocale("LC_CTYPE",".UTF-8"))
source(file.path(.current_root,"R/reanalysis_20261024.R"),encoding="UTF-8")
fl_current_entry <- function(args=character()) {{
  {note.strip()}
  pointer<-file.path(.current_root,"outputs/24_Two_generation_reanalysis/LATEST_RUN.txt")
  if(file.exists(pointer)&&!any(grepl("^--resume=",args)))args<-c(args,paste0("--resume=",readLines(pointer,n=1,warn=FALSE)))
  v3_main(args,stage="{stage}",root=.current_root)
}}
if(sys.nframe()==0L)fl_current_entry(commandArgs(TRUE))else if(isTRUE(getOption("fatty_liver.autorun",TRUE)))fl_current_entry()
'''
    (root/'scripts'/old.name).write_text(text,encoding='utf-8')
    manifest.append({'entry':old.name,'current_stage':stage,'original':str(old.relative_to(root))})
(root/'inventory/reanalysis_entry_migration_20261024.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print('Updated',len(manifest),'independent entries; old sources preserved.')
