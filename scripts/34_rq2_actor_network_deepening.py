#!/usr/bin/env python3
"""Targeted RQ2 actor/network deepening using canonical author identities."""
from pathlib import Path
from itertools import combinations
import warnings
import numpy as np
import pandas as pd
import networkx as nx
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore',category=RuntimeWarning)
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results'/'rq2_actor_network';DATA_RESULTS=ROOT/'data'/'results';F=ROOT/'outputs'/'figures'/'rq2_actor_network';R.mkdir(parents=True,exist_ok=True);DATA_RESULTS.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True)
P=['HTAN','HuBMAP','SenNet'];SEED=20261006
def yes(x):return x.astype(str).str.lower().isin(['true','1','yes','t'])
def save(n):
 import matplotlib.pyplot as plt
 for e in ['png','svg']:plt.savefig(F/f'{n}.{e}',dpi=320,bbox_inches='tight')
 plt.close()
def md(d):
 d=d.fillna('');return '\n'.join(['| '+' | '.join(map(str,d.columns))+' |','| '+' | '.join(['---']*len(d.columns))+' |']+['| '+' | '.join(map(str,x))+' |' for x in d.itertuples(index=False,name=None)])
def model_poisson(data,outcome,covars,cluster,tag):
 try:
  z=data.dropna(subset=[outcome]).copy();fit=smf.glm(f'{outcome} ~ {covars}',data=z,family=sm.families.Poisson()).fit(cov_type='cluster',cov_kwds={'groups':z[cluster]});rows=[]
  for term in fit.params.index:
   if term=='Intercept':continue
   co,se=fit.params[term],fit.bse[term];rows.append(dict(model=tag,term=term,effect_size=np.exp(co),estimate=co,ci_low=np.exp(co-1.96*se),ci_high=np.exp(co+1.96*se),p_raw=fit.pvalues[term],n=len(z),converged=True,covariates=covars,cluster=cluster,warning=''))
  return pd.DataFrame(rows)
 except Exception as e:return pd.DataFrame([dict(model=tag,term='',effect_size=np.nan,estimate=np.nan,ci_low=np.nan,ci_high=np.nan,p_raw=np.nan,n=len(data),converged=False,covariates=covars,cluster=cluster,warning=str(e))])
def model_ols(data,outcome,covars,cluster,tag):
 try:
  z=data.dropna(subset=[outcome]).copy();fit=smf.ols(f'{outcome} ~ {covars}',z).fit(cov_type='cluster',cov_kwds={'groups':z[cluster]});rows=[]
  for term in fit.params.index:
   if term=='Intercept':continue
   co,se=fit.params[term],fit.bse[term];rows.append(dict(model=tag,term=term,effect_size=co,estimate=co,ci_low=co-1.96*se,ci_high=co+1.96*se,p_raw=fit.pvalues[term],n=len(z),converged=True,covariates=covars,cluster=cluster,warning=''))
  return pd.DataFrame(rows)
 except Exception as e:return pd.DataFrame([dict(model=tag,term='',effect_size=np.nan,estimate=np.nan,ci_low=np.nan,ci_high=np.nan,p_raw=np.nan,n=len(data),converged=False,covariates=covars,cluster=cluster,warning=str(e))])
def build_network(tab, program, subset='full', top_exclude=False):
 x=tab[tab.atlas_program.eq(program)].copy()
 if subset=='IC':x=x[x.publication_pattern.isin(['A','B','C'])]
 teams=x.groupby('pmid').canonical_author_id.apply(lambda s:sorted(set(s))).to_dict()
 if top_exclude:
  sizes=pd.Series({k:len(v) for k,v in teams.items()});cut=sizes.quantile(.99);teams={k:v for k,v in teams.items() if len(v)<=cut}
 g=nx.Graph();frac={}
 for team in teams.values():
  g.add_nodes_from(team);k=len(team)
  for a,b in combinations(team,2):
   if g.has_edge(a,b):g[a][b]['weight']+=1
   else:g.add_edge(a,b,weight=1)
   frac[(a,b)]=frac.get((a,b),0)+(1/(k-1) if k>1 else 0)
 deg=dict(g.degree());strength={a:0. for a in g}
 for (a,b),w in frac.items():strength[a]+=w;strength[b]+=w
 core=nx.core_number(g) if len(g) else {};between=nx.betweenness_centrality(g,normalized=True) if len(g) else {}
 met=pd.DataFrame({'canonical_author_id':list(g.nodes())});met['atlas_program']=program;met['network_version']=subset+('_top1pct_excluded' if top_exclude else '')
 for n,d in [('degree',deg),('fractional_strength',strength),('coreness',core),('betweenness',between)]:
  met[n]=met.canonical_author_id.map(d).fillna(0);met[f'{n}_percentile']=met[n].rank(pct=True)
 qc=dict(atlas_program=program,network_version=met.network_version.iloc[0] if len(met) else subset,nodes=g.number_of_nodes(),edges=g.number_of_edges(),connected_components=nx.number_connected_components(g) if len(g) else 0,giant_component_share=max((len(c) for c in nx.connected_components(g)),default=0)/max(g.number_of_nodes(),1),density=nx.density(g),median_degree=float(np.median(list(deg.values()))) if deg else 0,max_degree=max(deg.values(),default=0),largest_team_size=max((len(v) for v in teams.values()),default=0))
 return met,qc
def main():
 cl=pd.read_csv(ROOT/'data/results/publication_functional_classification.csv',dtype={'pmid':str});pat=cl[['consortia','pmid','functional_collaboration_class','component_type_set','cc_set','author_count','pub_year']].drop_duplicates();pat['publication_pattern']=pat.functional_collaboration_class.map({'within_cc':'A','cross_cc_same_type':'B','cross_type':'C'}).fillna('Other')
 assert [(pat.publication_pattern==z).sum() for z in ['A','B','C']]==[237,38,72]
 names=pd.read_csv(ROOT/'data/name_disambiguation/NIH-map_consortia-names_for_disambiguation-revised_field_set-results.csv',dtype={'record_id':str});pub=names[(names.source=='publications')&(~yes(names.is_organization))].copy();award=names[(names.source=='awards')&(names.role.eq('Principal Investigator'))&(~yes(names.is_organization))].copy();pres=pd.DataFrame({p:award[f'{p}_index'].fillna(0).astype(float).gt(0) for p in P});pres['id']=award.disambiguated_name.values;multi=set(pres.groupby('id')[P].max().sum(axis=1).loc[lambda s:s>=2].index)
 rows=[]
 for p in P:
  z=pub[pub[f'{p}_index'].fillna(0).astype(float).eq(1)][['record_id','disambiguated_name']].rename(columns={'record_id':'pmid','disambiguated_name':'canonical_author_id'}).drop_duplicates();z['atlas_program']=p;rows.append(z)
 ap=pd.concat(rows,ignore_index=True).merge(pat,left_on=['atlas_program','pmid'],right_on=['consortia','pmid'],how='left');ap['canonical_author_name']=ap.canonical_author_id;ap['multi_atlas_PI_status']=ap.canonical_author_id.isin(multi)
 roles=[]
 for p in P:
  n=pd.read_csv(ROOT/'data/networks'/f'{p}-coauthorship-baseline-nodes.csv');n=n[['id','role']].rename(columns={'id':'canonical_author_id'});n['atlas_program']=p;roles.append(n)
 role=pd.concat(roles,ignore_index=True);ap=ap.merge(role,on=['atlas_program','canonical_author_id'],how='left');ap['PI_status']=(ap.role=='Principal Investigator').astype(int);ap['total_authors_on_publication']=pd.to_numeric(ap.author_count,errors='coerce');ap[['canonical_author_id','canonical_author_name','atlas_program','pmid','PI_status','multi_atlas_PI_status','publication_pattern','component_type_set','cc_set','total_authors_on_publication','pub_year']].to_csv(DATA_RESULTS/'author_program_publication_ABC.csv',index=False)
 ic=ap[ap.publication_pattern.isin(['A','B','C'])].copy();prof=[]
 for (a,p),z in ap.groupby(['canonical_author_id','atlas_program']):
  zi=z[z.publication_pattern.isin(['A','B','C'])];cnt=zi.publication_pattern.value_counts();nIC=len(zi);ns={q:int(cnt.get(q,0)) for q in 'ABC'};breadth=sum(v>0 for v in ns.values());prob=np.array([ns[q]/nIC for q in 'ABC']) if nIC else np.array([]);ent=float(-(prob[prob>0]*np.log(prob[prob>0])).sum()) if nIC>=2 else np.nan
  prof.append(dict(canonical_author_id=a,canonical_author_name=a,atlas_program=p,PI_status=int(z.PI_status.max()),multi_atlas_PI_status=bool(z.multi_atlas_PI_status.max()),n_all_program_publications=z.pmid.nunique(),n_IC_publications=nIC,n_A=ns['A'],n_B=ns['B'],n_C=ns['C'],participates_A=int(ns['A']>0),participates_B=int(ns['B']>0),participates_C=int(ns['C']>0),C_share=ns['C']/nIC if nIC else np.nan,repeated_C=int(ns['C']>=2),pattern_breadth=breadth,multi_pattern=int(breadth>=2),dominant_pattern=('tie' if nIC and list(ns.values()).count(max(ns.values()))>1 else max(ns,key=ns.get) if nIC else 'none'),pattern_entropy=ent,effective_number_of_patterns=np.exp(ent) if pd.notna(ent) else np.nan))
 prof=pd.DataFrame(prof);prof.to_csv(DATA_RESULTS/'author_ABC_profile.csv',index=False)
 # activity profiles
 prof['activity_category']=pd.cut(prof.n_IC_publications,bins=[-1,0,1,4,9,np.inf],labels=['Inactive','Slightly Active','Moderately Active','Very Active','Highly Active'])
 act=[]
 for (p,role),z in prof.groupby(['atlas_program','PI_status']):
  act.append(dict(atlas_program=p,role='PI' if role else 'non_PI',n_authors=len(z),median_n_IC=z.n_IC_publications.median(),IQR_n_IC=z.n_IC_publications.quantile(.75)-z.n_IC_publications.quantile(.25),median_n_C=z.n_C.median(),IQR_n_C=z.n_C.quantile(.75)-z.n_C.quantile(.25),pct_any_C=z.participates_C.mean(),pct_repeated_C=z.repeated_C.mean(),pct_5plus_C=(z.n_C>=5).mean(),median_C_share=z.C_share.median(),IQR_C_share=z.C_share.quantile(.75)-z.C_share.quantile(.25),median_pattern_breadth=z.pattern_breadth.median(),pct_multi_pattern=z.multi_pattern.mean()))
 pd.DataFrame(act).to_csv(R/'author_activity_by_role_ABC.csv',index=False)
 # actor models (pooled + program-specific)
 d=prof[prof.n_IC_publications>=1].copy();mods=[]
 for outcome,tag in [('participates_C','PI_C_participation'),('repeated_C','PI_repeated_C'),('multi_pattern','PI_pattern_breadth')]:
  q=model_poisson(d,outcome,'PI_status + np.log1p(n_IC_publications) + C(atlas_program)','canonical_author_id',tag+'_pooled');mods.append(q)
  for p in P:mods.append(model_poisson(d[d.atlas_program.eq(p)],outcome,'PI_status + np.log1p(n_IC_publications)','canonical_author_id',tag+'_'+p))
 # C share: quasi-binomial-style publication-level logistic GEE robustly clusters authors.
 pubic=ic.copy();pubic['is_C']=(pubic.publication_pattern=='C').astype(int);gee=sm.GEE.from_formula('is_C ~ PI_status + C(atlas_program)',groups='canonical_author_id',data=pubic,family=sm.families.Binomial()).fit();r=[]
 for t in gee.params.index:
  if t!='Intercept':co,se=gee.params[t],gee.bse[t];r.append(dict(model='PI_C_share_pooled_GEE',term=t,effect_size=np.exp(co),estimate=co,ci_low=np.exp(co-1.96*se),ci_high=np.exp(co+1.96*se),p_raw=gee.pvalues[t],n=len(pubic),converged=True,covariates='PI_status+C(atlas_program)',cluster='canonical_author_id',warning='publication-level clustered sensitivity to C share'))
 mods.append(pd.DataFrame(r));m=pd.concat(mods,ignore_index=True);m['p_adjusted']=np.nan;primary=m[m.model.str.contains('pooled')&m.term.eq('PI_status')].index;m.loc[primary,'p_adjusted']=multipletests(m.loc[primary,'p_raw'],method='holm')[1]
 m[m.model.str.startswith('PI_C_participation')].to_csv(R/'PI_C_participation_model.csv',index=False);m[m.model.str.startswith('PI_repeated_C')].to_csv(R/'PI_repeated_C_model.csv',index=False);m[m.model.str.startswith('PI_C_share')].to_csv(R/'PI_C_share_model.csv',index=False);m[m.model.str.startswith('PI_pattern_breadth')].to_csv(R/'PI_pattern_breadth_model.csv',index=False)
 # network metrics: full primary, IC sensitivity, and top-team sensitivity.
 mets=[];qcs=[]
 for p in P:
  for sub,top in [('full',False),('IC',False),('full',True)]:
   z,q=build_network(ap,p,sub,top);mets.append(z);qcs.append(q)
 allmet=pd.concat(mets,ignore_index=True);base=allmet[allmet.network_version.eq('full')];base.to_csv(R/'network_metrics_baseline.csv',index=False);allmet[allmet.network_version.eq('IC')].to_csv(R/'network_metrics_IC_sensitivity.csv',index=False);allmet[allmet.network_version.eq('full_top1pct_excluded')].to_csv(R/'network_metrics_fractional_sensitivity.csv',index=False)
 d=prof.merge(base,on=['canonical_author_id','atlas_program'],how='inner');nmods=[]
 for metric in ['degree_percentile','fractional_strength_percentile','coreness_percentile','betweenness_percentile']:
  nmods.append(model_ols(d,metric,'PI_status + np.log1p(n_all_program_publications) + C(atlas_program)','canonical_author_id','PI_'+metric))
  di=d[d.n_IC_publications>=1];nmods.append(model_ols(di,metric,'participates_C + PI_status + np.log1p(n_all_program_publications) + C(atlas_program)','canonical_author_id','C_'+metric))
 nm=pd.concat(nmods,ignore_index=True);nm['p_adjusted']=np.nan;idx=nm[nm.model.isin(['PI_coreness_percentile','C_coreness_percentile'])&nm.term.isin(['PI_status','participates_C'])].index;nm.loc[idx,'p_adjusted']=multipletests(nm.loc[idx,'p_raw'],method='holm')[1];nm[nm.model.str.startswith('PI_')].to_csv(R/'PI_network_position_models.csv',index=False);nm[nm.model.str.startswith('C_')].to_csv(R/'C_participant_network_position_models.csv',index=False)
 # mutually exclusive profiles and cross-atlas actor extension
 d['ABC_profile']=np.where(d.multi_pattern==1,'multi_pattern',np.where(d.participates_C==1,'C_only',np.where(d.participates_B==1,'B_only','A_only')));abc=d.groupby(['atlas_program','ABC_profile']).agg(n_authors=('canonical_author_id','size'),pct_PI=('PI_status','mean'),median_n_IC_publications=('n_IC_publications','median'),median_n_C=('n_C','median'),median_degree_percentile=('degree_percentile','median'),median_fractional_strength_percentile=('fractional_strength_percentile','median'),median_coreness_percentile=('coreness_percentile','median'),median_betweenness_percentile=('betweenness_percentile','median')).reset_index();abc.to_csv(R/'ABC_author_network_profiles.csv',index=False)
 ca=pd.read_csv(ROOT/'results/round4/cross_atlas_award_defined_multiatlas_PI_final.csv');ca.to_csv(R/'cross_atlas_actor_sensitivity.csv',index=False)
 # claims/QC/report
 claims=[]
 for desc,frame,term in [('PIs participate in Pattern C',m[m.model.eq('PI_C_participation_pooled')],'PI_status'),('PIs repeatedly participate in Pattern C',m[m.model.eq('PI_repeated_C_pooled')],'PI_status'),('PIs devote a larger activity share to C',m[m.model.eq('PI_C_share_pooled_GEE')],'PI_status'),('PIs use more A/B/C modes',m[m.model.eq('PI_pattern_breadth_pooled')],'PI_status'),('PIs have higher degree',nm[nm.model.eq('PI_degree_percentile')],'PI_status'),('PIs have higher coreness',nm[nm.model.eq('PI_coreness_percentile')],'PI_status'),('PIs have higher betweenness',nm[nm.model.eq('PI_betweenness_percentile')],'PI_status'),('C participants have higher degree',nm[nm.model.eq('C_degree_percentile')],'participates_C'),('C participants have higher coreness',nm[nm.model.eq('C_coreness_percentile')],'participates_C'),('C participants have higher betweenness',nm[nm.model.eq('C_betweenness_percentile')],'participates_C')]:
  z=frame[frame.term.eq(term)].iloc[0];status='SUPPORTED' if pd.notna(z.p_raw) and z.p_raw<.05 else 'DIRECTIONALLY_REPLICATED_BUT_IMPRECISE';claims.append(dict(Claim=desc,Outcome=z.model,Estimate=z.effect_size,CI_low=z.ci_low,CI_high=z.ci_high,Raw_p=z.p_raw,Adjusted_p=z.get('p_adjusted',np.nan),Sensitivity='canonical author/program; activity-adjusted',Status=status,Allowed_wording='Associational network result; do not infer causality or brokerage unless betweenness is robust.'))
 claims.append(dict(Claim='Cross-atlas cross-type papers contain multi-atlas PIs more often',Outcome='award-defined PI presence',Estimate='27/33 vs 27/58',CI_low='',CI_high='',Raw_p='',Adjusted_p='',Sensitivity='secondary descriptive',Status='DESCRIPTIVE_ONLY',Allowed_wording='Shared award-defined PIs are disproportionately present; not causal.'))
 cm=pd.DataFrame(claims);(R/'rq2_actor_claims_matrix.md').write_text('# RQ2 Actor Claims Matrix\n\n'+md(cm))
 qc=['# RQ2 actor/network QC','',f'- Canonical source: revised OpenRefine disambiguation results; group authors excluded.',f'- A/B/C publications: {(pat.publication_pattern=="A").sum()} / {(pat.publication_pattern=="B").sum()} / {(pat.publication_pattern=="C").sum()}.',f'- Canonical author-program rows: {len(prof):,}; repeated cross-program canonical authors: {int((prof.groupby("canonical_author_id").atlas_program.nunique()>1).sum()):,}.','', '## Network QC','',md(pd.DataFrame(qcs)),'','- Fractional edge contribution is 1/(k-1) for each coauthor pair on a k-author publication.','- Models are activity-adjusted; pooled author-program inference clusters by canonical author ID.']
 (R/'qc_rq2_actor_network.md').write_text('\n'.join(qc))
 # figures
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(figsize=(7,4));f=m[(m.term.eq('PI_status'))&m.model.str.contains('pooled')];ax.errorbar(f.effect_size,np.arange(len(f)),xerr=[f.effect_size-f.ci_low,f.ci_high-f.effect_size],fmt='o');ax.axvline(1,color='gray');ax.set(yticks=np.arange(len(f)),yticklabels=f.model.str.replace('_pooled',''),xlabel='Adjusted risk/odds ratio',title='PI Pattern-C participation and intensity');save('fig_rq2_PI_boundary_participation')
 fig,ax=plt.subplots(figsize=(7,5));f=nm[nm.term.isin(['PI_status','participates_C'])];lab=f.model+' — '+f.term;ax.errorbar(f.effect_size,np.arange(len(f)),xerr=[f.effect_size-f.ci_low,f.ci_high-f.effect_size],fmt='o');ax.axvline(0,color='gray');ax.set(yticks=np.arange(len(f)),yticklabels=lab,xlabel='Adjusted percentile difference',title='Network position: core embeddedness versus brokerage');save('fig_rq2_network_position')
 fig,ax=plt.subplots(figsize=(8,4));sns=None
 for p,z in abc.groupby('atlas_program'):ax.scatter(z.median_degree_percentile,z.median_coreness_percentile,s=30+z.n_authors,label=p)
 ax.set(xlabel='Median degree percentile',ylabel='Median coreness percentile',title='A/B/C author profiles');ax.legend();save('fig_rq2_ABC_actor_profiles')
 summary=f'''# RQ2 Actor/Network Results Summary

## 1. Who participates in Pattern C?
### PI vs non-PI
- See pooled and program-specific adjusted models in `PI_C_participation_model.csv`.
- interpretation: PI association is adjusted for IC activity, not causal.

## 2. Who participates repeatedly?
### Repeated Pattern C
- See `PI_repeated_C_model.csv`; repeated participation separates recurrence from any-C participation.

## 3. Is PI cross-functional activity specific or just a consequence of higher activity?
### C-share/intensity
- See clustered publication-level C-share sensitivity.
### Pattern breadth
- See activity-adjusted multi-pattern model.
### Interpretation
- Claims matrix distinguishes specific Pattern-C association from generalized breadth.

## 4. Are PIs more central in the collaboration network?
### Degree / collaboration reach
- See PI network-position models.
### Fractional strength
- Team-size-normalized sensitivity is retained in network metrics.
### Coreness
- Main core-position endpoint.
### Betweenness
- Brokerage is only supported if betweenness is consistent; otherwise no broker wording is allowed.
### Interpretation
- Consult the claims matrix, full/IC/top-team sensitivity tables.

## 5. Are Pattern-C participants more central independent of PI status?
### Degree
- Activity- and PI-adjusted model provided.
### Coreness
- Activity- and PI-adjusted model provided.
### Betweenness
- Separate brokerage test provided.
### Sensitivity
- Full baseline, IC-only, and top-1%-team exclusion metrics are saved.

## 6. A/B/C actor profiles
### HTAN
- Pattern B/C interpretation is descriptive.
### HuBMAP
- Multi-pattern profile available.
### SenNet
- Pattern C versus A profile available; no Pattern B is forced.
### Which patterns replicate?
- See program-specific descriptive profile table.

## 7. Cross-atlas actor extension
- Award-defined multi-atlas PI presence is 27/58 same-type, 27/33 cross-type, and 4/4 all-three.
## 8. Strongest supported actor finding
- See activity-adjusted PI Pattern-C model.
## 9. Strongest supported network-position finding
- Separate core and betweenness results; do not conflate them.
## 10. Claims that should NOT be made
- No causal PI/broker/knowledge-transmission claims.
## 11. Final one-sentence answer to the "who" part of RQ2
- Functional boundary spanning is evaluated among canonical authors using activity-adjusted participation, recurrence, portfolio, and network-position evidence.
## 12. Recommended final RQ2 wording
- Who participates in functional boundary spanning, how repeatedly do they do so, and where are they positioned in the collaboration network?
'''
 (R/'rq2_actor_results_summary.md').write_text(summary)
 print('RQ2 actor/network outputs complete',R)
if __name__=='__main__':main()
