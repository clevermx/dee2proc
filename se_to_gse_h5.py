%%time
import datetime
import bz2
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sps
from scipy.sparse import isspmatrix_csr,csr_matrix,csc_matrix,issparse,isspmatrix_csr
import math
cache_dir="D://datasets/DEE2/"
organism_name= 'athaliana'#'hsapiens' 'ecoli' mmusculus athaliana
geo_acc_name = 'GSE_accession' # may be "GSE_accession" or "Sample_name"
metadata_file = cache_dir+"metadata/"+ organism_name+"_meta_20200601.tsv"
data_file_name= cache_dir+"raw_files/star_counts/"+organism_name+"_se_20200530.tsv.bz2"
h5_gsm_name = cache_dir+"star_h5/"+organism_name+ "_star_matrix.h5"
gene_info_file_name =cache_dir+"GeneInfo/"+organism_name+"_GeneInfo.tsv"
now = datetime.datetime.now()
print("make gene index...")
gene_info=pd.read_csv(filepath_or_buffer=gene_info_file_name,sep="\t")
gene_ind= {gene:ind for ind,gene in enumerate(gene_info["GeneID"])}
n_genes=len(gene_ind)
print("gene quantity: "+str(n_genes))
print("read meta file...")
iter_meta = pd.read_csv(filepath_or_buffer=metadata_file,sep="\t", iterator=True, chunksize=1000)
meta_df = pd.concat([chunk for chunk in iter_meta])
true_srr=len(meta_df)
print("real SRR quantity: " +str(true_srr))
meta_df=meta_df[meta_df[geo_acc_name].str.startswith('GSM')]
n_srr= len(meta_df)
print("SRRs with gse: "+str(n_srr))
print("make SRR and gsm indexes...")
srr_to_gsm= dict(zip(meta_df["SRR_accession"],meta_df[geo_acc_name]))
srr_ind= {srr:ind for ind,srr in enumerate(srr_to_gsm.keys())}
gsm_ind= {gsm:ind for ind,gsm in enumerate(pd.unique(meta_df[geo_acc_name]))}
n_gsm= len(gsm_ind)
print("GSMs quantity:" )
gsm_from_srr_matrix= sps.lil_matrix((n_gsm,n_srr))
for srr_key in srr_to_gsm:
    gsm_from_srr_matrix[gsm_ind[srr_to_gsm[srr_key]],srr_ind[srr_key]]=1
print("start grouping...")
meta_df =meta_df.groupby(geo_acc_name).agg([lambda x: ';'.join(x)])
print("groupped")
meta_df=meta_df.loc[gsm_ind.keys()] 
print("reordered")
print("start h5...")
with h5py.File(h5_gsm_name, 'w') as h5_gse:
    print("create file and write meta...")
    meta_grp = h5_gse.create_group('meta')
    info_grp = h5_gse.create_group('info')
    data_grp= h5_gse.create_group('data')
    info_grp.create_dataset('version', data="dee2_gse_v1")
    info_grp.create_dataset('creation_date', data=now.strftime("%Y.%m.%d"))
    dt = h5py.special_dtype(vlen=np.str)
    meta_grp.create_dataset('genes',data=gene_ind.keys(),dtype=dt)      
    meta_grp.create_dataset('Sample_geo_accession',data=gsm_ind.keys(), dtype =dt)
    meta_grp.create_dataset('SRR_accession',data= np.array(meta_df["SRR_accession"], dtype=dt))
    meta_grp.create_dataset('Sample_instrument_model',data=np.array(meta_df["Model"], dtype=dt))
    meta_grp.create_dataset('Sample_last_update_date',data=np.array(meta_df["LoadDate"], dtype=dt))
    meta_grp.create_dataset('Sample_library_selection',data=np.array(meta_df["LibrarySelection"], dtype=dt))
    meta_grp.create_dataset('Sample_library_source',data=np.array(meta_df["LibrarySource"], dtype=dt))
    meta_grp.create_dataset('Sample_library_strategy',data=np.array(meta_df["LibraryStrategy"], dtype=dt))
    meta_grp.create_dataset('Sample_organism_ch1',data=np.array(meta_df["ScientificName"], dtype=dt))
    meta_grp.create_dataset('Sample_series_id',data=np.array(meta_df["Samplename"], dtype=dt))
    meta_grp.create_dataset('Sample_status',data=np.array(meta_df["Consent"], dtype=dt))
    meta_grp.create_dataset('Sample_submission_date',data=np.array(meta_df["ReleaseDate"], dtype=dt))
    meta_grp.create_dataset('Sample_quality',data=np.array(meta_df["QC_summary"], dtype=dt))
    meta_grp.create_dataset('Sample_title',data=np.array(meta_df["BioSample"], dtype=dt))
    meta_grp.create_dataset('Sample_type',data=np.full(len(meta_df["BioSample"]),"SRA", dtype=dt))
    exp_data=data_grp.create_dataset("expression", (n_gsm,n_genes),dtype= 'i4')#, compression="gzip", compression_opts=9)
    srr_per_time=1000
    iter_data = pd.read_csv(filepath_or_buffer=data_file_name,sep="\t", iterator=True, chunksize=srr_per_time*n_genes,names=["srr","gene","se"])
    proc_srr= 0
    gsm_from_srr_matrix=csr_matrix(gsm_from_srr_matrix)    
    print("start process expression...")
    for chunk in iter_data:
        global_ids =[]
        local_ids=[]
        local_srr_ind= {srr:ind for ind,srr in enumerate(pd.unique(chunk['srr']))}
        srr_count=len(local_srr_ind)
        proc_srr=proc_srr+srr_count
        for cur_srr in local_srr_ind:
            try:
                srr_id= srr_ind[cur_srr]
                global_ids.append(srr_id)
                local_ids.append(local_srr_ind[cur_srr])
            except(KeyError):
                continue
        if local_ids:
            raw_matrix= chunk["se"].values.astype(int).reshape(srr_count,n_genes) #13.6 милисек
            A=gsm_from_srr_matrix[:,global_ids]  #0.135 милисек
            gsm_mask= [np.sum(x)>0 for x in A] #4 милисек
            A = A.tocsr()  #очень быстро            
            B=csr_matrix(raw_matrix[local_ids,]) #9 милисек
            A = A.dot(B) #6 милисек
            exp_data[gsm_mask== True,0:n_genes]+= A[gsm_mask== True,0:n_genes] #425 милисек
        if proc_srr%500 ==0:
            print(proc_srr/true_srr)
    h5_gse.flush()
print(1)

  