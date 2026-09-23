"""Pre-specified sensitivity grid for the ex vivo / clinical concordance."""
import pandas as pd, numpy as np, itertools
from scipy import stats
warn=[]
pr=pd.read_csv('pairs_recist.csv'); pf=pd.read_csv('pairs_pfs.csv')
qc=pd.read_csv('qc_by_patient.csv').set_index('study_id').qc_verdict.to_dict()
pr['qc']=pr.study_id.map(lambda s: qc.get(s,'UNKNOWN'))
CLASS_TRI='class-level (trametinib ex vivo / cobimetinib clinical)'

rows=[]
for cutoff in (30,32,50):
  for excl_exh in (True,False):
    for qc_gate in (True,False):
      for excl_class in (True,False):
        for partial_resp in (True,False):
          d=pr.copy()
          if excl_exh:   d=d[~d.flag]
          if qc_gate:    d=d[d.qc!='FAIL']
          if excl_class: d=d[d.class_match!=CLASS_TRI]
          if len(d)<3: 
              rows.append(dict(cutoff=cutoff,exclude_tissue_exhausted=excl_exh,qc_gate=qc_gate,
                               exclude_class_level=excl_class,partial_counts_as_response=partial_resp,
                               n_pairs=len(d),n_patients=d.study_id.nunique(),
                               TP=np.nan,FP=np.nan,FN=np.nan,TN=np.nan,accuracy=np.nan,
                               fisher_p=np.nan,spearman_rho=np.nan)); continue
          # ex vivo call at this cutoff
          d['sens_call']=d.pct_of_control < (100-cutoff)
          # clinical responder definition
          resp_set={'CR','PR'} | ({'Partial'} if partial_resp else set())
          d['resp']=d.response.isin(resp_set)
          TP=int((d.sens_call&d.resp).sum()); FP=int((d.sens_call&~d.resp).sum())
          FN=int((~d.sens_call&d.resp).sum()); TN=int((~d.sens_call&~d.resp).sum())
          acc=(TP+TN)/len(d)*100
          try: fp_=stats.fisher_exact([[TP,FP],[FN,TN]])[1]
          except Exception: fp_=np.nan
          rho=stats.spearmanr(d.pct_of_control,d.recist_rank)[0] if d.pct_of_control.nunique()>1 else np.nan
          rows.append(dict(cutoff=cutoff,exclude_tissue_exhausted=excl_exh,qc_gate=qc_gate,
                           exclude_class_level=excl_class,partial_counts_as_response=partial_resp,
                           n_pairs=len(d),n_patients=d.study_id.nunique(),
                           TP=TP,FP=FP,FN=FN,TN=TN,accuracy=round(acc,1),
                           fisher_p=round(fp_,3) if fp_==fp_ else np.nan,
                           spearman_rho=round(rho,3) if rho==rho else np.nan))
g=pd.DataFrame(rows)
g.to_csv('sensitivity_grid.csv',index=False)
pd.set_option('display.width',250)
print("=== SENSITIVITY GRID: ex vivo call vs RECIST responder ===")
print(g.to_string(index=False))
print(f"\nscenarios: {len(g)}")
sub=g.dropna(subset=['accuracy'])
print(f"n_pairs range {sub.n_pairs.min()}-{sub.n_pairs.max()}, patients {sub.n_patients.min()}-{sub.n_patients.max()}")
print(f"accuracy range {sub.accuracy.min():.0f}-{sub.accuracy.max():.0f}%")
print(f"Spearman rho range {sub.spearman_rho.min():+.2f} to {sub.spearman_rho.max():+.2f}")
print(f"scenarios with Fisher p < 0.05: {(sub.fisher_p<0.05).sum()} of {len(sub)}")
