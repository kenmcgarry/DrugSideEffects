library(xtable)

common<-c(
"Fatigue","Asthenia","Hypotension","Confusion","Postural hypotension","Cough","Haemorrhage",
"Dyspepsia","Anxiety","Diarrhoea","Abdominal pain","Malaise","Dry mouth","Gastroenteritis",
"Syncope","Constipation","Blurred vision","Flatulence","Somnolence","Dizziness","Headache",
"Paresthesia","Anorexia","Vomiting","Insomnia","Chest pain","Rhinitis","Back pain","Anxiety","Weight loss","Tinnitus","","")

uncommon <- c("Urinary retention","Arrhythmia","Hallucinations","Tremor","Agitation","Increased salivation", 
              "Fever","Epistaxis","Urinary tract infection","Hypertension","Anemia","Incontinence", 
              "Gastroenteritis","Apathy","Supraventricular tachycardia","Pneumonia","leg cramps","Melena",
              "Arthritis","Infection","Hypertonia","Vertigo","AV block","Bronchitis","Ischemia","Heart failure", 
              "Nervousness","libido increased","Convulsions","Dysphagia","Rash","Ataxia","Hyperglycemia")
              
  
rare <- c(
  "Nausea","renal failure","Myalgia","Gastritis","Rectal hemorrhage","Bradycardia","Atrial fibrillation","Hypokalemia",  
  "Nocturia","Hypersensitivity","Alkaline phosphatase increased","Oedema","Fibrillation","Hypokinesia",  
  "Urinary System Disorders","Psychiatric Disorders","Transient ischemic attack","Aggressive","Dehydration","Thrombocytopenia",
  "Peripheral edema","Arthralgia","Cramps","Gastrointestinal haemorrhage","Urinary incontinence","Cerebrovascular accident", 
  "Aphasia","Hematuria","Cataract","Tachycardia","Urinary frequency","Diverticulitis","Infarction")  


df <- data.frame(common = common,uncommon = uncommon,rare=rare)
df <- print.xtable(xtable(df),include.rownames=FALSE)
print.xtable(df,label='side97',caption = 'The 97 combined side-effects for the three Alzhiemers drugs organised into common, uncommon and rare cases')



