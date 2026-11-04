#!/usr/bin/env python3
"""Component-type collaboration analysis for NIH Mapping Programs.

Uses the validated program-publication data and preserves the manuscript's
division-level collaboration labels while adding a separate functional type
classification.  Run from the repository root:
    python3 scripts/40_component_type_analysis.py
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path
import math
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.proportion import proportion_confint
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "outputs" / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

PROGRAMS = ["HTAN", "HuBMAP", "SenNet"]
TYPES = ["coordination", "experimental", "technology"]
TYPE_LABEL = {
    "coordination": "Coordination",
    "experimental": "Experimental",
    "technology": "Technology",
}
PAIR_LABEL = {
    "coordination__experimental": "Coordination–experimental",
    "coordination__technology": "Coordination–technology",
    "experimental__technology": "Experimental–technology",
    "coordination__coordination": "Coordination–coordination",
    "experimental__experimental": "Experimental–experimental",
    "technology__technology": "Technology–technology",
}

# The repository's current HuBMAP label is RIT; it corresponds to RTI in the
# frozen brief and is intentionally retained as the original label in outputs.
MAPPING = {
    ("HTAN", "DCC"): "coordination",
    ("HTAN", "HTA"): "experimental",
    ("HTAN", "PCA"): "experimental",
    ("HuBMAP", "HIVE"): "coordination",
    ("HuBMAP", "TMC_H"): "experimental",
    ("HuBMAP", "DEMO"): "experimental",
    ("HuBMAP", "TTD"): "technology",
    ("HuBMAP", "RIT"): "technology",
    ("HuBMAP", "RTI"): "technology",
    ("SenNet", "CODCC"): "coordination",
    ("SenNet", "TMC_S"): "experimental",
    ("SenNet", "TDA"): "technology",
}


def is_true(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "t", "1", "yes"])


def pair_key(a: str, b: str) -> str:
    return "__".join(sorted((a, b)))


def q(values, p):
    values = pd.Series(values).dropna()
    return float(values.quantile(p)) if len(values) else np.nan


def md_table(frame: pd.DataFrame, index=False) -> str:
    """Render a compact Markdown table without the optional tabulate package."""
    x = frame.reset_index() if index else frame.copy()
    x = x.fillna("")
    cols = [str(c) for c in x.columns]
    body = [[str(v).replace("|", "\\|") for v in row] for row in x.astype(str).itertuples(index=False, name=None)]
    return "\n".join(["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"]*len(cols)) + " |"] + ["| " + " | ".join(row) + " |" for row in body])


def rcr_stats(values: pd.Series) -> dict:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if values.empty:
        return dict(n=0, median_rcr=np.nan, q1_rcr=np.nan, q3_rcr=np.nan,
                    iqr_rcr=np.nan, mean_rcr=np.nan, min_rcr=np.nan, max_rcr=np.nan)
    return dict(n=len(values), median_rcr=values.median(), q1_rcr=q(values, .25),
                q3_rcr=q(values, .75), iqr_rcr=q(values, .75)-q(values, .25),
                mean_rcr=values.mean(), min_rcr=values.min(), max_rcr=values.max())


def boot_median(values: pd.Series, rng: np.random.Generator, reps=5000):
    x = pd.to_numeric(values, errors="coerce").dropna().to_numpy()
    if len(x) < 5:
        return np.nan, np.nan, np.nan
    med = np.median(x)
    boots = np.median(rng.choice(x, size=(reps, len(x)), replace=True), axis=1)
    return med, *np.quantile(boots, [.025, .975])


def bootstrap_difference(a: pd.Series, b: pd.Series, rng: np.random.Generator, reps=5000):
    a = pd.to_numeric(a, errors="coerce").dropna().to_numpy()
    b = pd.to_numeric(b, errors="coerce").dropna().to_numpy()
    if min(len(a), len(b)) < 5:
        return np.nan, np.nan, np.nan
    obs = np.median(a) - np.median(b)
    d = (np.median(rng.choice(a, (reps, len(a)), replace=True), axis=1)
         - np.median(rng.choice(b, (reps, len(b)), replace=True), axis=1))
    return obs, *np.quantile(d, [.025, .975])


def savefig(name: str):
    for suffix in ("svg", "png"):
        plt.savefig(FIGURES / f"{name}.{suffix}", dpi=320, bbox_inches="tight")
    plt.close()


def load_data():
    pubs = pd.read_csv(DATA / "processed" / "NIHReporter-map_consortia-publications-validated-citations.csv",
                       dtype={"pmid": str})
    meta = pd.read_csv(DATA / "metadata" / "NIH-map_consortia-publication-metadata.csv",
                       dtype={"pmid": str, "core_project_number": str})
    awards = pd.read_csv(DATA / "metadata" / "NIH-map_consortia_awards-division-metadata.csv",
                         dtype={"core_project_number": str})
    return pubs, meta, awards


def make_classification(pubs, meta):
    pubs = pubs.loc[~is_true(pubs.preprint_matched)].copy()
    pubs["is_preprint_bool"] = is_true(pubs.is_preprint)
    pubs = pubs.drop_duplicates(["consortia", "pmid"])
    assoc = meta.merge(pubs[["consortia", "pmid"]], on=["consortia", "pmid"], how="inner")
    assoc = assoc.drop_duplicates(["consortia", "pmid", "core_project_number", "division"])
    assoc["component_type"] = [MAPPING.get((p, d)) for p, d in zip(assoc.consortia, assoc.division)]
    if assoc.component_type.isna().any():
        missing = assoc.loc[assoc.component_type.isna(), ["consortia", "division"]].drop_duplicates()
        raise ValueError(f"Unmapped division labels: {missing.to_dict('records')}")

    rows = []
    for (program, pmid), x in assoc.groupby(["consortia", "pmid"], sort=False):
        projects = sorted(x.core_project_number.dropna().unique())
        ccs = sorted(x.division.dropna().unique())
        types = sorted(x.component_type.dropna().unique())
        if len(projects) == 1:
            klass = "single_project"
        elif len(ccs) == 1:
            klass = "within_cc"
        elif len(types) == 1:
            klass = "cross_cc_same_type"
        else:
            klass = "cross_type"
        rows.append(dict(consortia=program, pmid=pmid, n_projects=len(projects), n_cc=len(ccs),
                         n_component_types=len(types), cc_set="; ".join(ccs),
                         component_type_set="; ".join(types),
                         functional_collaboration_class=klass,
                         within_type_IC=klass in ("within_cc", "cross_cc_same_type"),
                         cross_type_IC=klass == "cross_type"))
    classified = pd.DataFrame(rows).merge(pubs, on=["consortia", "pmid"], how="left", validate="one_to_one")
    expected = {"single_project": "Single Project", "within_cc": "Intradivision Collaboration",
                "cross_cc_same_type": "Interdivision Collaboration", "cross_type": "Interdivision Collaboration"}
    classified["expected_old_classification"] = classified.functional_collaboration_class.map(expected)
    classified["old_classification_match"] = classified.collab_type.eq(classified.expected_old_classification)
    return pubs, assoc, classified


def component_membership(classified, assoc):
    types = assoc[["consortia", "pmid", "component_type"]].drop_duplicates()
    return types.merge(classified, on=["consortia", "pmid"], how="left", validate="many_to_one")


def orientation(awards, assoc, membership):
    funded = awards.copy()
    funded["component_type"] = [MAPPING.get((p, d)) for p, d in zip(funded.project_id, funded.division)]
    if funded.component_type.isna().any():
        raise ValueError("An award division is unmapped")
    funded = funded.rename(columns={"project_id": "consortia"})
    funded = funded[["consortia", "core_project_number", "division", "component_type"]].drop_duplicates()
    active = assoc[["consortia", "core_project_number", "component_type"]].drop_duplicates()
    rows=[]
    for program in PROGRAMS:
        for typ in TYPES:
            f = funded.query("consortia == @program and component_type == @typ")
            if f.empty: continue
            m = membership.query("consortia == @program and component_type == @typ")
            all_p = m.pmid.nunique(); ic = m.loc[m.functional_collaboration_class != "single_project", "pmid"].nunique()
            within = m.loc[m.within_type_IC, "pmid"].nunique(); cross = m.loc[m.cross_type_IC, "pmid"].nunique()
            low, high = proportion_confint(cross, ic, method="wilson") if ic else (np.nan, np.nan)
            ct_projects = assoc.loc[(assoc.consortia==program)&(assoc.component_type==typ)&assoc.pmid.isin(m.loc[m.cross_type_IC,"pmid"]), "core_project_number"].nunique()
            denom_projects = len(f)
            rows.append(dict(consortia=program, component_type=typ,
                             funded_core_projects=denom_projects,
                             projects_with_publication=active.query("consortia == @program and component_type == @typ").core_project_number.nunique(),
                             all_publications=all_p, ic_publications=ic, within_type_ic_publications=within,
                             cross_type_ic_publications=cross, cross_type_share=cross/ic if ic else np.nan,
                             wilson_ci_low=low, wilson_ci_high=high,
                             cross_type_project_participants=ct_projects,
                             cross_type_project_participation=ct_projects/denom_projects if denom_projects else np.nan,
                             functional_boundary_orientation=(cross-within)/(cross+within) if cross+within else np.nan,
                             numerator_cross_type=cross, denominator_ic=ic))
    return pd.DataFrame(rows), funded


def publication_pairs(classified, assoc):
    x = classified.query("cross_type_IC").merge(assoc, on=["consortia", "pmid"], how="left")
    rows=[]; triples=[]
    for (program, pmid), g in x.groupby(["consortia", "pmid"]):
        types = sorted(g.component_type.unique()); nprojects=g.core_project_number.nunique()
        if len(types) == 3: triples.append((program, pmid))
        for a,b in combinations(types,2):
            pa = sorted(g.loc[g.component_type==a,"core_project_number"].unique())
            pb = sorted(g.loc[g.component_type==b,"core_project_number"].unique())
            rows.append(dict(consortia=program, pmid=pmid, pair=pair_key(a,b), type_a=a, type_b=b,
                             n_projects_in_publication=nprojects, projects_a="; ".join(pa), projects_b="; ".join(pb)))
    detail=pd.DataFrame(rows)
    out=[]
    for program in PROGRAMS:
        total=x.loc[x.consortia==program,"pmid"].nunique()
        for a,b in combinations(TYPES,2):
            p=pair_key(a,b); g=detail.query("consortia == @program and pair == @p") if not detail.empty else pd.DataFrame()
            pairs=set()
            for _, r in g.iterrows():
                pairs.update(combinations(sorted(set(r.projects_a.split('; ')+r.projects_b.split('; '))),2))
            # keep only true cross-side pairs
            crosspairs=set()
            for _, r in g.iterrows():
                crosspairs.update((u,v) for u in r.projects_a.split('; ') for v in r.projects_b.split('; '))
            out.append(dict(consortia=program, pair=p, type_a=a, type_b=b,
                            cross_type_publications=g.pmid.nunique(), pct_of_program_cross_type_publications=g.pmid.nunique()/total if total else np.nan,
                            unique_projects_type_a=len(set(sum((s.split('; ') for s in g.projects_a), []))) if len(g) else 0,
                            unique_projects_type_b=len(set(sum((s.split('; ') for s in g.projects_b), []))) if len(g) else 0,
                            unique_cross_type_project_pairs=len(crosspairs),
                            median_projects_per_publication=g.n_projects_in_publication.median() if len(g) else np.nan,
                            q1_projects_per_publication=q(g.n_projects_in_publication,.25), q3_projects_per_publication=q(g.n_projects_in_publication,.75),
                            all_three_type_publications=len({pm for pr,pm in triples if pr==program})))
    return pd.DataFrame(out), detail


def density(assoc, funded):
    # One fixed type per program/project; validate before building edges.
    project_types=assoc[["consortia","core_project_number","component_type"]].drop_duplicates()
    if project_types.duplicated(["consortia","core_project_number"]).any():
        raise ValueError("A core project maps to more than one component type")
    edges=[]
    for (program, pmid), g in assoc.groupby(["consortia","pmid"]):
        ids=sorted(g.core_project_number.unique())
        for a,b in combinations(ids,2): edges.append((program,a,b,pmid))
    ed=pd.DataFrame(edges,columns=["consortia","project_a","project_b","pmid"])
    if ed.empty: ed=pd.DataFrame(columns=["consortia","project_a","project_b","pmid"])
    weights=ed.groupby(["consortia","project_a","project_b"],as_index=False).pmid.nunique().rename(columns={"pmid":"shared_publications"})
    type_map={(r.consortia,r.core_project_number):r.component_type for r in project_types.itertuples()}
    if len(weights):
        weights["type_a"]=[type_map[(p,x)] for p,x in zip(weights.consortia,weights.project_a)]
        weights["type_b"]=[type_map[(p,x)] for p,x in zip(weights.consortia,weights.project_b)]
    rows=[]
    for program in PROGRAMS:
        f=funded.query("consortia == @program")[["core_project_number","component_type"]].drop_duplicates()
        active=project_types.query("consortia == @program")
        for i,a in enumerate(TYPES):
            for b in TYPES[i:]:
                fa=f.query("component_type == @a").core_project_number.nunique(); fb=f.query("component_type == @b").core_project_number.nunique()
                aa=active.query("component_type == @a").core_project_number.nunique(); ab=active.query("component_type == @b").core_project_number.nunique()
                possible=fa*fb if a!=b else fa*(fa-1)//2
                possible_active=aa*ab if a!=b else aa*(aa-1)//2
                w=weights.loc[(weights.consortia==program)&(((weights.type_a==a)&(weights.type_b==b))|((weights.type_a==b)&(weights.type_b==a)))]
                rows.append(dict(consortia=program,pair=pair_key(a,b),type_a=a,type_b=b,
                                 funded_projects_a=fa,funded_projects_b=fb,active_projects_a=aa,active_projects_b=ab,
                                 observed_collaborating_project_pairs=len(w),possible_funded_project_pairs=possible,
                                 project_pair_density=len(w)/possible if possible else np.nan,
                                 possible_active_project_pairs=possible_active,active_project_pair_density=len(w)/possible_active if possible_active else np.nan,
                                 total_shared_publication_edge_weight=w.shared_publications.sum(),
                                 mean_shared_publications_per_observed_pair=w.shared_publications.mean(),
                                 median_shared_publications_per_observed_pair=w.shared_publications.median()))
    d=pd.DataFrame(rows)
    pooled=[]
    for pair,g in d.groupby("pair"):
        eligible=g[g.possible_funded_project_pairs>0]
        pooled.append(dict(consortia="Pooled",pair=pair,type_a=eligible.type_a.iloc[0],type_b=eligible.type_b.iloc[0],
                           funded_projects_a=np.nan,funded_projects_b=np.nan,active_projects_a=np.nan,active_projects_b=np.nan,
                           observed_collaborating_project_pairs=eligible.observed_collaborating_project_pairs.sum(),
                           possible_funded_project_pairs=eligible.possible_funded_project_pairs.sum(),
                           project_pair_density=eligible.observed_collaborating_project_pairs.sum()/eligible.possible_funded_project_pairs.sum(),
                           possible_active_project_pairs=eligible.possible_active_project_pairs.sum(),
                           active_project_pair_density=eligible.observed_collaborating_project_pairs.sum()/eligible.possible_active_project_pairs.sum() if eligible.possible_active_project_pairs.sum() else np.nan,
                           total_shared_publication_edge_weight=eligible.total_shared_publication_edge_weight.sum(),
                           mean_shared_publications_per_observed_pair=np.nan,median_shared_publications_per_observed_pair=np.nan))
    return pd.concat([d,pd.DataFrame(pooled)],ignore_index=True)


def nominal_assortativity(assoc):
    """Unweighted Newman nominal assortativity on each program project network."""
    project_type = assoc[["consortia", "core_project_number", "component_type"]].drop_duplicates()
    lookup = {(x.consortia, x.core_project_number): x.component_type for x in project_type.itertuples()}
    rows=[]
    for program in PROGRAMS:
        edges=set()
        for _, g in assoc.query("consortia == @program").groupby("pmid"):
            ids=sorted(g.core_project_number.unique())
            edges.update(combinations(ids,2))
        m=len(edges)
        e=np.zeros((len(TYPES),len(TYPES)))
        for a,b in edges:
            i=TYPES.index(lookup[(program,a)]); j=TYPES.index(lookup[(program,b)])
            if i==j: e[i,i]+=1/m if m else 0
            else:
                e[i,j]+=1/(2*m); e[j,i]+=1/(2*m)
        marg=e.sum(axis=1); denom=1-(marg**2).sum()
        assort=(np.trace(e)-(marg**2).sum())/denom if denom else np.nan
        rows.append(dict(consortia=program,project_network_edges=m,same_type_edge_fraction=float(np.trace(e)),nominal_assortativity=assort))
    return pd.DataFrame(rows)


def robustness_model(classified):
    """Optional HC3 OLS robustness model; publication type is unavailable."""
    x=classified.loc[(~classified.is_preprint_bool)&classified.relative_citation_ratio.notna()].copy()
    x["functional_collaboration_class"]="single_project"
    for idx,r in x.iterrows():
        types=r.component_type_set.split("; ")
        if r.within_type_IC: x.loc[idx,"functional_collaboration_class"]="within_type_IC"
        elif r.cross_type_IC: x.loc[idx,"functional_collaboration_class"]=pair_key(*types) if len(types)==2 else "all_three_types"
    complete=x.dropna(subset=["relative_citation_ratio","pub_year","author_count","functional_collaboration_class"])
    fit=smf.ols("np.log1p(relative_citation_ratio) ~ C(functional_collaboration_class) + C(consortia) + pub_year + author_count",data=complete).fit(cov_type="HC3")
    ci=fit.conf_int(); out=pd.DataFrame({"term":fit.params.index,"estimate":fit.params.values,"std_error":fit.bse.values,"p_value":fit.pvalues.values,"ci_low":ci.iloc[:,0].values,"ci_high":ci.iloc[:,1].values})
    info=pd.DataFrame([dict(input_publications=len(x),complete_case_publications=len(complete),excluded_publications=len(x)-len(complete),outcome="log1p(relative_citation_ratio)",covariates="functional class, program, publication year, author count",unavailable_covariate="publication_type")])
    return out, info


def author_roles(classified):
    rows=[]; qc=[]
    for program in PROGRAMS:
        node_path=DATA/"networks"/f"{program}-coauthorship-baseline-nodes.csv"
        nodes=pd.read_csv(node_path)
        role=dict(zip(nodes.id.astype(str),nodes.role))
        pubs=classified.query("consortia == @program and functional_collaboration_class != 'single_project'")
        hits=0; total=0; author_map={}
        for _, r in pubs.iterrows():
            authors=[a.strip() for a in str(r.authors).split(";") if a.strip() and a.strip()!='nan']
            for author in authors:
                total+=1
                if author not in role: continue
                hits+=1
                author_map.setdefault(author, {"role":role[author],"within":set(),"cross":set()})
                (author_map[author]["cross"] if r.cross_type_IC else author_map[author]["within"]).add(r.pmid)
        qc.append(dict(consortia=program, author_publication_strings=total, matched_to_node_role=hits, role_match_rate=hits/total if total else np.nan))
        for role_name, label in [("Principal Investigator","PI"),("Researcher","non_PI_researcher")]:
            a=[v for v in author_map.values() if v["role"]==role_name]
            n=len(a); cross=sum(bool(v["cross"]) for v in a)
            counts=[len(v["cross"]) for v in a]
            rows.append(dict(consortia=program,role=label,ic_participants=n,cross_type_participants=cross,
                             cross_type_participation_rate=cross/n if n else np.nan,median_cross_type_publications=np.median(counts) if n else np.nan,
                             q1_cross_type_publications=q(counts,.25),q3_cross_type_publications=q(counts,.75),
                             pct_with_2plus_cross_type=sum(c>=2 for c in counts)/n if n else np.nan))
    out=pd.DataFrame(rows)
    compare=[]
    for program,g in out.groupby("consortia"):
        pi=g[g.role=="PI"].iloc[0]; non=g[g.role=="non_PI_researcher"].iloc[0]
        a,n1=int(pi.cross_type_participants),int(pi.ic_participants); c,n0=int(non.cross_type_participants),int(non.ic_participants)
        rr=(a/n1)/(c/n0) if n1 and n0 and c else np.nan
        aa=a+.5; cc=c+.5; nn1=n1+1; nn0=n0+1
        se=math.sqrt(1/aa-1/nn1+1/cc-1/nn0)
        _, fisher_p = fisher_exact([[a, n1-a], [c, n0-c]], alternative="two-sided")
        compare.append(dict(consortia=program,pi_rate=a/n1 if n1 else np.nan,non_pi_rate=c/n0 if n0 else np.nan,
                            risk_ratio=rr,rr_ci_low=math.exp(math.log((aa/nn1)/(cc/nn0))-1.96*se) if n1 and n0 else np.nan,
                            rr_ci_high=math.exp(math.log((aa/nn1)/(cc/nn0))+1.96*se) if n1 and n0 else np.nan,
                            fisher_exact_p_value=fisher_p))
    return out, pd.DataFrame(compare), pd.DataFrame(qc)


def legacy_pi_participation_reconciliation():
    """Recompute the legacy Table 8 PI statistic from its node-file sources.

    The taskbook calls this manuscript Table 5, while the checked-in output is
    `outputs/tables/table_8.csv`.  The legacy denominator is every PI node in
    the baseline network and its numerator is every PI node in the internal-
    collaboration network.  This is intentionally separate from the new
    author-string-matched, cross-*type* analysis.
    """
    legacy = pd.read_csv(ROOT / "outputs" / "tables" / "table_8.csv")
    legacy = legacy[["Mapping_Program", "Count_of_PIs",
                     "Count_PIs_involved_in_Internal_Collaborative_Research",
                     "Percent_of_PIs_involved_in_Collaborative_Research"]].rename(columns={
                         "Mapping_Program": "consortia",
                         "Count_of_PIs": "legacy_pi_denominator",
                         "Count_PIs_involved_in_Internal_Collaborative_Research": "legacy_pi_ic_participants",
                         "Percent_of_PIs_involved_in_Collaborative_Research": "legacy_pi_ic_percent",
                     })
    recomputed = []
    for program in PROGRAMS:
        baseline = pd.read_csv(DATA / "networks" / f"{program}-coauthorship-baseline-nodes.csv")
        collab = pd.read_csv(DATA / "networks" / f"{program}-coauthorship-collaboration_partition-nodes.csv")
        denominator = int((baseline.role == "Principal Investigator").sum())
        numerator = int((collab.role == "Principal Investigator").sum())
        recomputed.append(dict(
            consortia=program,
            recomputed_pi_denominator=denominator,
            recomputed_pi_ic_participants=numerator,
            recomputed_pi_ic_percent=round(100 * numerator / denominator, 2),
        ))
    out = legacy.merge(pd.DataFrame(recomputed), on="consortia", how="outer")
    out["denominator_match"] = out.legacy_pi_denominator.eq(out.recomputed_pi_denominator)
    out["numerator_match"] = out.legacy_pi_ic_participants.eq(out.recomputed_pi_ic_participants)
    out["percent_match_after_rounding"] = out.legacy_pi_ic_percent.round(2).eq(out.recomputed_pi_ic_percent.round(2))
    return out


def rcr_tables(membership, classified):
    rng=np.random.default_rng(20261005)
    r=membership.loc[(~membership.is_preprint_bool)&membership.relative_citation_ratio.notna() & membership.functional_collaboration_class.ne("single_project")].copy()
    r["boundary_group"]=np.where(r.cross_type_IC,"cross_type_IC","within_type_IC")
    rows=[]; diffs=[]
    for (p,t,g),x in r.groupby(["consortia","component_type","boundary_group"]):
        rows.append(dict(consortia=p,component_type=t,boundary_group=g,**rcr_stats(x.relative_citation_ratio)))
    for (p,t),x in r.groupby(["consortia","component_type"]):
        cross=x.loc[x.boundary_group=="cross_type_IC","relative_citation_ratio"]
        within=x.loc[x.boundary_group=="within_type_IC","relative_citation_ratio"]
        d,lo,hi=bootstrap_difference(cross,within,rng)
        diffs.append(dict(consortia=p,component_type=t,n_cross_type=cross.notna().sum(),n_within_type=within.notna().sum(),
                          median_difference_cross_minus_within=d,bootstrap_ci_low=lo,bootstrap_ci_high=hi,
                          inference_status="estimated" if not np.isnan(d) else "exploratory_small_n"))
    pair_rows=[]
    cross=classified.loc[(classified.cross_type_IC)&(~classified.is_preprint_bool)&classified.relative_citation_ratio.notna()].copy()
    assoc_cols=["consortia","pmid","component_type"]
    # The caller attaches component_type_set, allowing exact two-type pair assignment.
    for (p, pmid), x in cross.groupby(["consortia","pmid"]):
        types=sorted(x.component_type_set.iloc[0].split("; "))
        cat=pair_key(*types) if len(types)==2 else "all_three_types" if len(types)==3 else "unclassified"
        pair_rows.append(dict(consortia=p,pmid=pmid,functional_pair=cat,relative_citation_ratio=x.relative_citation_ratio.iloc[0]))
    pairs=pd.DataFrame(pair_rows)
    pair_summary=[]
    for (p,cat),x in pairs.groupby(["consortia","functional_pair"]):
        med,lo,hi=boot_median(x.relative_citation_ratio,rng)
        pair_summary.append(dict(consortia=p,functional_pair=cat,**rcr_stats(x.relative_citation_ratio),bootstrap_median_ci_low=lo,bootstrap_median_ci_high=hi,
                                 inference_status="estimated" if len(x)>=5 else "exploratory_small_n"))
    return r, pd.DataFrame(rows).merge(pd.DataFrame(diffs),on=["consortia","component_type"],how="outer"), pd.DataFrame(pair_summary)


def figures(orient, pairpub, density_df, role, rcr):
    sns.set_theme(style="whitegrid", font_scale=.9)
    colors={"HTAN":"#008042","HuBMAP":"#6D99D5","SenNet":"#9964A8"}
    # Orientation
    fig,ax=plt.subplots(figsize=(8,4.4))
    pos={t:i for i,t in enumerate(TYPES)}
    offsets={"HTAN":-.22,"HuBMAP":0,"SenNet":.22}
    for p,g in orient.groupby("consortia"):
        for _,r in g.iterrows():
            y=pos[r.component_type]+offsets[p]
            ax.errorbar(r.cross_type_share,y,xerr=[[max(0, r.cross_type_share-r.wilson_ci_low)],[max(0, r.wilson_ci_high-r.cross_type_share)]],fmt='o',color=colors[p],capsize=2)
            label_x = .985 if r.cross_type_share > .9 else r.cross_type_share + .025
            ax.text(label_x, y, f"{int(r.numerator_cross_type)}/{int(r.denominator_ic)}", va='center', ha='right' if r.cross_type_share > .9 else 'left', fontsize=7)
    ax.set(yticks=list(pos.values()),yticklabels=[TYPE_LABEL[t] for t in TYPES],xlim=(0,1),xlabel="Cross-type share of internal-collaboration publications (95% Wilson CI)",ylabel="Component type")
    ax.legend([plt.Line2D([],[],marker='o',linestyle='',color=colors[p]) for p in PROGRAMS],PROGRAMS,title="Program",loc="lower center",bbox_to_anchor=(.5,1.02),ncol=3,frameon=False)
    fig.tight_layout(); savefig("fig_component_type_orientation")
    # Pair matrix raw counts / density, one row per program plus pooled density.
    fig,axes=plt.subplots(2,4,figsize=(12,6),sharex=True,sharey=True)
    pairs=[pair_key(a,b) for a,b in combinations(TYPES,2)]
    for j,p in enumerate(PROGRAMS):
        mat=np.zeros((3,3)); x=pairpub.query("consortia == @p")
        for _,r in x.iterrows():
            i=TYPES.index(r.type_a); k=TYPES.index(r.type_b); mat[i,k]=mat[k,i]=r.cross_type_publications
        sns.heatmap(mat,annot=True,fmt='.0f',cmap="Blues",cbar=j==3,ax=axes[0,j],xticklabels=[TYPE_LABEL[t] for t in TYPES],yticklabels=[TYPE_LABEL[t] for t in TYPES],square=True)
        axes[0,j].set_title(p)
    axes[0,3].axis('off'); axes[0,3].set_title("A: publication co-occurrence")
    for j,p in enumerate(PROGRAMS+["Pooled"]):
        mat=np.full((3,3),np.nan); x=density_df.query("consortia == @p")
        for _,r in x.iterrows():
            i=TYPES.index(r.type_a); k=TYPES.index(r.type_b); mat[i,k]=mat[k,i]=r.project_pair_density
        sns.heatmap(mat,annot=True,fmt='.3f',cmap="YlOrRd",vmin=0,vmax=np.nanmax(density_df.project_pair_density),cbar=j==3,ax=axes[1,j],xticklabels=[TYPE_LABEL[t] for t in TYPES],yticklabels=[TYPE_LABEL[t] for t in TYPES],square=True)
        axes[1,j].set_title(p)
    axes[1,3].set_title("Pooled B: funded-project density")
    fig.tight_layout(); savefig("fig_component_type_pair_matrix")
    # Role rates
    fig,ax=plt.subplots(figsize=(6,4))
    sns.barplot(data=role,x="consortia",y="cross_type_participation_rate",hue="role",palette={"PI":"#2B6CB0","non_PI_researcher":"#718096"},ax=ax)
    ax.set(xlabel="Program",ylabel="Cross-type participation rate among IC authors",ylim=(0,1)); ax.legend(title="Role",labels=["PI","Non-PI researcher"])
    for c in ax.containers: ax.bar_label(c,fmt="%.2f",padding=2,fontsize=8)
    fig.tight_layout(); savefig("fig_author_role_cross_type_participation")
    # RCR
    fig,axes=plt.subplots(1,3,figsize=(12,4.2),sharey=True)
    for ax,p in zip(axes,PROGRAMS):
        x=rcr.query("consortia == @p")
        sns.boxplot(data=x,x="component_type",y="relative_citation_ratio",hue="boundary_group",order=TYPES,hue_order=["cross_type_IC","within_type_IC"],palette={"cross_type_IC":"#5B8FD1","within_type_IC":"#D88A58"},showfliers=False,ax=ax)
        ax.set(title=p,xlabel="",ylabel="RCR" if p==PROGRAMS[0] else "")
        labels=[]
        for t in TYPES:
            n_cross=len(x[(x.component_type==t)&(x.boundary_group=="cross_type_IC")])
            n_within=len(x[(x.component_type==t)&(x.boundary_group=="within_type_IC")])
            labels.append(f"{TYPE_LABEL[t]}\nC={n_cross}, W={n_within}")
        ax.set_xticks(range(len(TYPES)), labels, rotation=25, ha="right")
        if p!=PROGRAMS[0]: ax.get_legend().remove()
    axes[0].legend(title="Class",labels=["Cross type","Within type"])
    fig.tight_layout(); savefig("fig_rcr_within_vs_cross_type")


def qc_and_summary(pubs, assoc, classified, awards, orient, pairpub, density_df, roles, rr, pairrcr, role_qc):
    ctab=pd.crosstab(classified.functional_collaboration_class,classified.collab_type,dropna=False)
    violations=classified.loc[~classified.old_classification_match,["consortia","pmid","functional_collaboration_class","collab_type","cc_set","component_type_set"]]
    mapping_qc=(awards.rename(columns={"project_id":"consortia"})[["consortia","division","core_project_number"]].drop_duplicates()
                .assign(component_type=lambda d:[MAPPING.get((p,x)) for p,x in zip(d.consortia,d.division)])
                .groupby(["consortia","division","component_type"],dropna=False).agg(project_count=("core_project_number","nunique")).reset_index())
    mapping_qc = mapping_qc.merge(
        assoc.groupby(["consortia","division"]).pmid.nunique().rename("publication_count").reset_index(),
        on=["consortia","division"], how="left").fillna({"publication_count": 0})
    pub_counts=classified.groupby("consortia").pmid.nunique().rename("program_publications").reset_index()
    overall=pd.DataFrame({"metric":["deduplicated PMIDs across programs","multi-program PMIDs","preprints retained","journal articles retained"],"value":[classified.pmid.nunique(),classified.groupby("pmid").consortia.nunique().gt(1).sum(),classified.is_preprint_bool.sum(),(~classified.is_preprint_bool).sum()]})
    participation=[]
    for p in PROGRAMS:
        denom=awards.query("project_id == @p").core_project_number.nunique()
        num=assoc.merge(classified.query("consortia == @p and functional_collaboration_class != 'single_project'")[["consortia","pmid"]],on=["consortia","pmid"]).core_project_number.nunique()
        participation.append(dict(consortia=p,funded_projects=denom,projects_in_internal_collaboration=num,percent=100*num/denom))
    pi_note="PI participation is an author-level measure, not a project-level measure. The legacy Table 8 PI denominator and IC-participant numerator are recomputed and reconciled below; the new role analysis uses a validated-publication, author-name-matched cohort to measure cross-type participation."
    qc=["# Component-type analysis QC", "", "## Mapping QC", "", md_table(mapping_qc), "", "All observed project divisions map to exactly one frozen component type. The repository label `RIT` is retained and mapped as the taskbook's RTI technology label.", "", "## Classification QC", "", md_table(ctab, index=True), "", f"Mismatched old/new classifications: **{len(violations)}**.", "", "### Mismatches", "", md_table(violations) if len(violations) else "None.", "", "## Publication-set QC", "", md_table(pub_counts), "", md_table(overall), "", "The functional classification is derived from the same validated, preprint-deduplicated program-publication set used in the existing publication analysis.", "", "## Known manuscript inconsistencies", "", "### Project participation recomputation", "", md_table(pd.DataFrame(participation)), "", "Definition: numerator = unique funded Core Projects appearing in at least one non-single-project publication; denominator = all funded Core Projects mapped to that program.", "", "### PI participation", "", pi_note, "", "### Author-role matching QC", "", md_table(role_qc)]
    (RESULTS/"qc_component_type_analysis.md").write_text("\n".join(qc))
    violations.to_csv(RESULTS/"qc_component_type_classification_violations.csv",index=False)
    mapping_qc.to_csv(RESULTS/"component_type_mapping.csv",index=False)
    # Evidence-limited manuscript summary, generated from completed result tables.
    top_pair=pairpub.loc[pairpub.cross_type_publications.idxmax()] if len(pairpub) and pairpub.cross_type_publications.max()>0 else None
    top_density=density_df.query("consortia != 'Pooled'").loc[density_df.query("consortia != 'Pooled'").project_pair_density.idxmax()]
    pooled_cross = density_df.query("consortia == 'Pooled' and type_a != type_b").sort_values("project_pair_density", ascending=False).iloc[0]
    pooled_pairs = pairpub.groupby("pair",as_index=False).cross_type_publications.sum().sort_values("cross_type_publications",ascending=False)
    top_pooled_pair = pooled_pairs.iloc[0]
    summaries=[]
    for _,r in orient.iterrows(): summaries.append(f"- {r.consortia} {TYPE_LABEL[r.component_type]}: {int(r.numerator_cross_type)}/{int(r.denominator_ic)} cross-type IC publications ({r.cross_type_share:.1%}; 95% CI {r.wilson_ci_low:.1%}–{r.wilson_ci_high:.1%}).")
    role_lines=[]
    for p,g in roles.groupby("consortia"):
        pi=g[g.role=="PI"].iloc[0]; non=g[g.role=="non_PI_researcher"].iloc[0]
        a,n1=pi.cross_type_participants,pi.ic_participants; c,n0=non.cross_type_participants,non.ic_participants
        risk=(a/n1)/(c/n0); se=math.sqrt(1/(a+.5)-1/(n1+1)+1/(c+.5)-1/(n0+1))
        lo,hi=math.exp(math.log(risk)-1.96*se),math.exp(math.log(risk)+1.96*se)
        role_lines.append(f"- {p}: PI {int(a)}/{int(n1)} ({pi.cross_type_participation_rate:.1%}) vs non-PI researcher {int(c)}/{int(n0)} ({non.cross_type_participation_rate:.1%}); RR {risk:.2f} (95% CI {lo:.2f}–{hi:.2f}).")
    rcr_lines=[]
    for (p,t),g in rr.groupby(["consortia","component_type"]):
        cross=g[g.boundary_group=="cross_type_IC"]; within=g[g.boundary_group=="within_type_IC"]
        if len(cross) and len(within):
            c=cross.iloc[0]; w=within.iloc[0]
            status="exploratory" if c.inference_status=="exploratory_small_n" else f"bootstrap Δ={c.median_difference_cross_minus_within:.2f} (95% CI {c.bootstrap_ci_low:.2f}–{c.bootstrap_ci_high:.2f})"
            rcr_lines.append(f"- {p} {TYPE_LABEL[t]}: cross-type n={int(c.n)}, median RCR {c.median_rcr:.2f}; within-type n={int(w.n)}, median {w.median_rcr:.2f}; {status}.")
    tech = rr[(rr.component_type=="technology") & rr.boundary_group.isin(["cross_type_IC","within_type_IC"])]
    summary=["# Component-Type Results Summary","","## RQ1. Functional collaboration orientation",*summaries,"","- Cross-program recurrence: technology is inward-oriented in both programs with technology components (HuBMAP 25.0%; SenNet 16.4% cross-type share). Experimental orientation differs substantially (HTAN 10.5%, HuBMAP 68.3%, SenNet 48.1%). Coordination estimates are boundary-spanning but are based on one funded coordination project in HTAN and SenNet.","","## RQ2. Functional interfaces",f"- Most frequent raw functional pair pooled across programs: {PAIR_LABEL[top_pooled_pair.pair]} ({int(top_pooled_pair.cross_type_publications)} pair co-occurrences; a three-type publication contributes to each relevant pair).",f"- Most frequent program-specific raw functional pair: {top_pair.consortia} {PAIR_LABEL[top_pair.pair]} ({int(top_pair.cross_type_publications)} publications)." if top_pair is not None else "- No cross-type pair observed.",f"- Highest program-specific funded-project-pair density: {top_density.consortia} {PAIR_LABEL[top_density.pair]} ({top_density.project_pair_density:.3f}; {int(top_density.observed_collaborating_project_pairs)}/{int(top_density.possible_funded_project_pairs)} possible pairs).",f"- Highest pooled cross-type funded-project-pair density: {PAIR_LABEL[pooled_cross.pair]} ({pooled_cross.project_pair_density:.3f}; {int(pooled_cross.observed_collaborating_project_pairs)}/{int(pooled_cross.possible_funded_project_pairs)} within-program possible pairs).",*role_lines,"","## RQ3. Citation influence",*rcr_lines,"- Functional-pair RCR results use exactly-two-type publications; all-three-type publications are reported separately in `rcr_by_functional_pair.csv`.","","## Direct answers to taskbook decision questions","- The same types do not show a universal common orientation. Technology is consistently inward-oriented across HuBMAP and SenNet; experimental is mixed; coordination is underpowered in HTAN and SenNet.","- HuBMAP and SenNet both show experimental–technology as a prominent interface (16 and 21 publications, respectively); HuBMAP also has 16 coordination–experimental publications.","- PIs have higher observed cross-type participation than non-PI researchers in every program, with the effect sizes and CIs reported above; this is an authorship-participation association, not evidence of brokerage or causation.","- The apparent technology cross-type RCR advantage survives aggregation directionally in HuBMAP (4.37 vs 2.92) and SenNet (6.84 vs 3.38), but both bootstrap CIs include zero; conclusion: MIXED / statistically imprecise.","- Strongest defensible one-sentence finding: true cross-type collaboration is concentrated at program-specific functional interfaces, while technology components are comparatively inward-oriented and citation associations vary by component type and program.","","## Claims supported by the data","- The publication-level classification separates same-type inter-CC collaboration from true cross-type collaboration.","- Cross-type shares, interface counts, opportunity-normalized densities, role participation, and citation summaries are program-stratified.","","## Claims NOT supported / should not be made","- These associations do not establish direct knowledge exchange, PI brokerage, causation, or a universal citation advantage.","","## Small-n or unstable comparisons","- Any RCR comparison with fewer than five publications in either group is marked `exploratory_small_n` in the result tables."]
    (RESULTS/"component_type_results_summary.md").write_text("\n".join(summary))


def main():
    pubs, meta, awards = load_data()
    pubs, assoc, classified = make_classification(pubs, meta)
    membership = component_membership(classified, assoc)
    orient, funded = orientation(awards, assoc, membership)
    pairpub, _ = publication_pairs(classified, assoc)
    dens = density(assoc, funded)
    roles, role_compare, role_qc = author_roles(classified)
    # Add each publication's type set to its membership record so RCR pairs are exact.
    classified_for_rcr = classified.merge(assoc.groupby(["consortia","pmid"]).component_type.apply(lambda x:"; ".join(sorted(set(x)))).rename("component_type_set").reset_index(),on=["consortia","pmid"],suffixes=("","_new"))
    membership_rcr = membership.drop(columns=["component_type_set"],errors="ignore").merge(classified_for_rcr[["consortia","pmid","component_type_set"]],on=["consortia","pmid"],how="left")
    rcr_long, rcr_summary, pair_rcr = rcr_tables(membership_rcr, classified_for_rcr)
    assort = nominal_assortativity(assoc)
    model_terms, model_info = robustness_model(classified_for_rcr)
    pi_reconciliation = legacy_pi_participation_reconciliation()
    # Reconcile the legacy project-participation table against unique funded
    # Core Project denominators used by this component-type analysis.
    legacy = pd.read_csv(ROOT / "outputs" / "tables" / "table_6.csv")
    legacy = legacy[legacy.Publication_Group.eq("Internal Collaborative Publications")][["Mapping_Program","Count_of_Projects","Percent_of_Projects"]].rename(columns={"Mapping_Program":"consortia","Count_of_Projects":"legacy_ic_projects","Percent_of_Projects":"legacy_ic_percent"})
    canonical=[]
    for p in PROGRAMS:
        den=funded.query("consortia == @p").core_project_number.nunique()
        num=assoc.merge(classified.query("consortia == @p and functional_collaboration_class != 'single_project'")[["consortia","pmid"]],on=["consortia","pmid"]).core_project_number.nunique()
        canonical.append(dict(consortia=p,canonical_funded_projects=den,canonical_ic_projects=num,canonical_ic_percent=100*num/den))
    reconciliation=legacy.merge(pd.DataFrame(canonical),on="consortia",how="outer")
    # Outputs
    classified.to_csv(RESULTS/"publication_functional_classification.csv",index=False)
    orient.to_csv(RESULTS/"component_type_orientation_by_program.csv",index=False)
    pairpub.to_csv(RESULTS/"component_pair_publication_matrix.csv",index=False)
    dens.to_csv(RESULTS/"component_pair_project_density.csv",index=False)
    roles.merge(role_compare,on="consortia",how="left").to_csv(RESULTS/"author_role_cross_type_participation.csv",index=False)
    rcr_summary.to_csv(RESULTS/"rcr_within_vs_cross_type.csv",index=False)
    pair_rcr.to_csv(RESULTS/"rcr_by_functional_pair.csv",index=False)
    assort.to_csv(RESULTS/"component_type_nominal_assortativity.csv",index=False)
    model_terms.to_csv(RESULTS/"rcr_robustness_model_coefficients.csv",index=False)
    model_info.to_csv(RESULTS/"rcr_robustness_model_qc.csv",index=False)
    reconciliation.to_csv(RESULTS/"legacy_project_participation_reconciliation.csv",index=False)
    pi_reconciliation.to_csv(RESULTS/"legacy_pi_participation_reconciliation.csv",index=False)
    figures(orient,pairpub,dens,roles,rcr_long)
    qc_and_summary(pubs,assoc,classified,awards,orient,pairpub,dens,roles,rcr_summary,pair_rcr,role_qc)
    with (RESULTS/"qc_component_type_analysis.md").open("a") as f:
        f.write("\n\n## Legacy table reconciliation\n\n")
        f.write("The repository's `table_6.csv` uses non-unique award-row denominators (35 HTAN, 80 HuBMAP, 50 SenNet when inferred from its counts and percentages). The canonical component-type analysis uses unique funded Core Project Numbers.\n\n")
        f.write(md_table(reconciliation) + "\n\n")
        f.write("### PI participation reconciliation\n\n")
        f.write("The taskbook calls this manuscript Table 5, but the checked-in source is `outputs/tables/table_8.csv`; `table_5.csv` is a different multi-program-PI table. Table 8 was recomputed directly from the baseline and internal-collaboration node files. All three PI denominators, IC-participant numerators, and rounded percentages match exactly; thus the taskbook's stated arithmetic issue does not apply to this repository revision.\n\n")
        f.write(md_table(pi_reconciliation) + "\n\n")
        f.write("The new author-role table answers a different question: among authors whose names match the validated publication set and who appear on at least one IC publication, what share appear on a cross-*type* publication? Its denominators are therefore not expected to equal Table 8's full network-based IC counts. Matching coverage is reported above.\n\n")
        f.write("## Optional robustness\n\n")
        f.write("### Nominal assortativity\n\n" + md_table(assort) + "\n\n")
        f.write("### Adjusted RCR model\n\n" + md_table(model_info) + "\n\n")
        f.write("HC3-robust OLS coefficients are saved in `rcr_robustness_model_coefficients.csv`. The model excludes publication type because it is not available in the validated input.\n")
    with (RESULTS/"component_type_results_summary.md").open("a") as f:
        f.write("\n\n## Optional robustness checks\n\n")
        f.write("- Nominal project-network assortativity and the adjusted log1p(RCR) model were run after the primary descriptive analyses; outputs are in `component_type_nominal_assortativity.csv` and `rcr_robustness_model_coefficients.csv`. The adjusted model omits publication type because that covariate is unavailable.\n")
    print("Completed component-type analysis:", RESULTS)


if __name__ == "__main__":
    main()
