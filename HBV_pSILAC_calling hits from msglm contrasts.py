#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd #version 2.2.2
import numpy as np #version 1.26.4
import matplotlib.pyplot as plt #version 3.9.0
import matplotlib.colors as mcolors
import seaborn as sns # version 0.13.2
import pickle


# In[3]:


contrasts=pd.read_csv("object_contrasts.csv") #load the msglm output

contrasts_unfiltered=contrasts.copy()
# In[5]:


#basic filtering

def filter_count(df, nrows, nproteins):
    nrows.append(len(df))
    nproteins.append(df["object_id"].nunique())

    print(
        f"total rows: {len(df)},\n"
        f"total proteins: {df['object_id'].nunique()}"
    )
    return nrows, nproteins

nrows=[]
nproteins=[]
nrows, nproteins = filter_count(contrasts, nrows, nproteins)
contrasts=contrasts[contrasts["std_type"]=="median"]
nrows,nproteins = filter_count(contrasts, nrows, nproteins)
contrasts=contrasts[contrasts["is_contaminant"]==False]
nrows,nproteins = filter_count(contrasts, nrows, nproteins)
contrasts=contrasts[contrasts["is_reverse"]==False]
nrows,nproteins = filter_count(contrasts, nrows, nproteins)
contrasts=contrasts[contrasts["timepoint_lhs"]==contrasts["timepoint_rhs"]] #only matching timepoints
nrows,nproteins = filter_count(contrasts, nrows, nproteins)
contrasts=contrasts[contrasts["treatment_lhs"]!=contrasts["treatment_rhs"]] #only differing treatments
nrows,nproteins = filter_count(contrasts, nrows, nproteins)
contrasts=contrasts[contrasts["treatment_rhs"]!="HBV"] #exclude HBV vs HBVdX contrast; only compare vs mock
nrows,nproteins = filter_count(contrasts, nrows, nproteins)
contrasts=contrasts[(contrasts["nmsruns_quanted_lhs_min"]>=3)|(contrasts["nmsruns_quanted_rhs_min"]>=3)]  #3 replicates quantified on either side for a valid contrast
nrows,nproteins = filter_count(contrasts, nrows, nproteins)


# In[ ]:


#acccess data coverage for H and L channel by timepoint
df = contrasts.copy()#[~((contrasts["timepoint_lhs"] == 0) & (contrasts["mstag"] == "H"))]

df_counts = df.groupby(
    ["treatment_lhs", "mstag", "object_id"]
)["timepoint_lhs"].nunique().reset_index()

#export numbers from mstag timepoint count to excel:
summary = (
    df_counts
    .groupby(["treatment_lhs", "mstag", "timepoint_lhs"])
    .agg(protein_count=("object_id", "nunique"))
    .reset_index()
)
summary["percent"] = summary.groupby(
    ["treatment_lhs", "mstag"]
)["protein_count"].transform(lambda x: 100 * x / x.sum())
table = summary.pivot_table(
    index=["treatment_lhs", "mstag"],
    columns="timepoint_lhs",
    values=["protein_count", "percent"],
    fill_value=0
)
table.to_excel("timepoint_distribution_by_virus_mstag.xlsx")


# In[ ]:


#acess identified proteins by timepoint
# proteins per timepoint
df = contrasts[~((contrasts["timepoint_lhs"] == 0) & (contrasts["mstag"] == "H"))]
tp_count=df.groupby(["treatment_lhs","timepoint_lhs"])["object_id"].nunique().reset_index()
tp_count.to_excel("Protein_count_per_timepoint.xlsx")


# In[11]:


df = contrasts[~((contrasts["timepoint_lhs"] == 0) & (contrasts["mstag"] == "H"))]
contrasts["major protein_acs"]=contrasts['majority_protein_acs'].str.split(";").str[0]

#pickle a dictionary mapping object_label to UniProt to use in other scripts
contrasts[["object_label","major protein_acs"]]
ObjectLabel_to_UniProt=dict(zip(contrasts["object_label"], contrasts["major protein_acs"]))
with open("ObjectLabel_to_UniProt.pkl", "wb") as f:  # 'wb' = write binary
    pickle.dump(ObjectLabel_to_UniProt, f)
    
#create a dataframe where I have the H and L channel in same row
contrasts["contrast w/o tag"]=contrasts["contrast"].str.split(":").str[1]
dfH=contrasts[contrasts["mstag"]=="H"]
dfL=contrasts[contrasts["mstag"]=="L"]
DF=dfH.merge(dfL, on=["object_id","contrast w/o tag"], how="outer", suffixes=("_H","_L"))


# In[12]:


#identifying synthesis and degradation

#strong degredation candidates
DF["deg"]=((DF["prob_nonneg_L"] < 0.005) #if L is likely negative
           & (DF["median_log2_L"] < -0.5) #if L is decreasing strongly
           & (DF["p_value_L"] < 0.005)
           & (DF["median_log2_H"] > -0.1)) # if L decrease is not due to reduced synthesis (kind of assuming that taking synthesis at t is good reflection of t-1)

#strong stabilisation candidates
DF["stab"]=((DF["prob_nonpos_L"] < 0.005) #if L is definetly positive
           & (DF["median_log2_L"] > 0.5) #if L is strongy up
           & (DF["p_value_L"] < 0.005)
           & (DF["median_log2_H"] < 0.2)) #and this cannot be explained through increased synthesis

#strong synthesis up candidates
DF["syn_up"]=((DF["prob_nonpos_H"] < 0.005) #if H is definetly positive
              & (DF["median_log2_H"] > 0.5)
              & (DF["p_value_H"] < 0.005)) 

#strong synthesis down candidates
DF["syn_down"]=((DF["prob_nonneg_H"] < 0.005) #if H is definetly negative
              & (DF["p_value_H"] < 0.005)
              & (DF["median_log2_H"] < -0.5)) #if H shows a strong decrease


# In[13]:


#clean up dataframe to be able to call a timepoint, object_label or treatment_lhs (HBV or HBVdX; treatment_rhs always mock) without any of the nan from merging
DF["object_label"]=DF["object_label_L"].fillna(DF["object_label_H"])
DF["timepoint"]=DF["timepoint_lhs_L"].fillna(DF["timepoint_lhs_H"])
DF["treatment_lhs"]=DF["treatment_lhs_L"].fillna(DF["treatment_lhs_H"])


# In[14]:


for i in ["deg","stab","syn_up","syn_down"]:
    df=DF[DF[i]==True]
    num_tp=len(df)
    prots=df["object_id"].unique()
    print(f"for {i} there are {len(prots)} proteins found and {num_tp} significant contrasts (single time point)")


# In[16]:


#only accept changes if they occured in at least 2 consecutive time points

#get columns of either deg, stab, syn_up, syn_down == True (work on one at a time)
#group DF by object_id (not object_label, seems to be not unique), mstag, and treatment lhs (treatment_rhs is alsways mock)
# then collect timepoints into one column and check if [1,2],[2,3],[3,4],[4,6],[6,8] or [8,10] is in there.
# if not, use object_id, Treatment_lhs, mstag to set the column of change ==False 


allowed_pairs = {(1,2), (2,3), (3,4), (4,6), (6,8), (8,10)}
tp_collapsed_all = [] 

for change in ["deg","stab","syn_up","syn_down"]:
    print(f"proteins which are hit for {change} before selection for two consecutive timepoints: {DF[DF[change]==True]["object_id"].nunique()}")
    tp_collapsed = (
        DF[DF[change]]
        .groupby(["object_id", "treatment_lhs"])["timepoint"]
        .apply(lambda x: sorted(set(x)))
        .reset_index(name="timepoints")
    )
    
    tp_collapsed["two_consec_tp"] = tp_collapsed["timepoints"].apply(
        lambda tp: any(
            (a, b) in allowed_pairs
            for a, b in zip(tp[:-1], tp[1:])
        )
    )
    
    tp_collapsed["regulation"] = change
    
    print(f"proteins found to check:{tp_collapsed["object_id"].nunique()}")
    print(f"proteins pass check: {tp_collapsed[tp_collapsed["two_consec_tp"]==True]["object_id"].nunique()}")
    
    keys = ["object_id", "treatment_lhs"]
    
    mapper = (
        tp_collapsed
        .set_index(keys)["two_consec_tp"]
    )
    tp_collapsed_all.append(tp_collapsed)

    keys = ["object_id", "treatment_lhs"]

    mapper = tp_collapsed.set_index(keys)["two_consec_tp"]

    # map to DF (temporary Series)
    two_consec = (
        DF.set_index(keys).index.map(mapper).fillna(False)
    )

    # only invalidate existing hits
    DF.loc[DF[change] & ~two_consec, change] = False

    print(f"proteins after: {DF.loc[DF[change], 'object_id'].nunique()}")

# FINAL combined table
tp_collapsed = pd.concat(tp_collapsed_all, ignore_index=True)


# In[18]:


#save data of all proteins to explore
DF.to_excel("HBV_pSILAC_DF_all_proteins.xlsx")


# In[19]:


#filter by two-consecutive timepoint criteria to get stable regulations=="hits"
tp_collapsed[tp_collapsed["timepoints"].apply(lambda x: len(x))>1]
tp_collapsed["timepoints"] = tp_collapsed["timepoints"].apply(lambda x: sorted(x))
object_id_to_label_mapper=dict(zip(DF["object_id"],DF["object_label"]))
tp_collapsed["object_label"]=tp_collapsed["object_id"].map(object_id_to_label_mapper)
Hits_summary=tp_collapsed[tp_collapsed["two_consec_tp"]]
Hits_summary.to_excel("HBV_pSILAC_hits_summary.xlsx")

