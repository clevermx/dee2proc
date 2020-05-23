import bz2
import h5py
import numpy as np
import pandas as pd
import datetime
organism_name= 'ecoli'
cache_dir="D://datasets/DEE2/"
h5_file_name = cache_dir+"/srr_files/"+organism_name+ "_SRR_matrix.h5"
meta_file_name = cache_dir+"/metadata/"+organism_name+"_metadata.tsv"
data_file_name = cache_dir +"/raw_files/"+organism_name +"_se.tsv.bz2"
with h5py.File(h5_file_name, 'w') as h5_file:
    meta_grp = h5_file.create_group('meta')
    info_grp = h5_file.create_group('info')
    data_grp= h5_file.create_group('data')
    info_grp.create_dataset('version', data="dee2_v1")
    info_grp.create_dataset('creation_date', data="2020-02-23")
    meta_df=pd.read_csv(filepath_or_buffer=meta_file_name,sep="\t")
    dt = h5py.special_dtype(vlen=np.str)
    meta_grp.create_dataset('SRR_accession',data= np.array(meta_df["SRR_accession"], dtype=dt))
    meta_grp.create_dataset('Sample_geo_accession',data=np.array(meta_df["GSE_accession"], dtype=dt))
    meta_grp.create_dataset('Sample_instrument_model',data=np.array(meta_df["Model"], dtype=dt))
    meta_grp.create_dataset('Sample_last_update_date',data=np.array(meta_df["LoadDate"], dtype=dt))
    meta_grp.create_dataset('Sample_library_selection',data=np.array(meta_df["LibrarySelection"], dtype=dt))
    meta_grp.create_dataset('Sample_library_source',data=np.array(meta_df["LibrarySource"], dtype=dt))
    meta_grp.create_dataset('Sample_library_strategy',data=np.array(meta_df["LibraryStrategy"], dtype=dt))
    meta_grp.create_dataset('Sample_organism_ch1',data=np.array(meta_df["ScientificName"], dtype=dt))
    meta_grp.create_dataset('Sample_series_id',data=np.array(meta_df["Sample_name"], dtype=dt))
    meta_grp.create_dataset('Sample_status',data=np.array(meta_df["Consent"], dtype=dt))
    meta_grp.create_dataset('Sample_submission_date',data=np.array(meta_df["ReleaseDate"], dtype=dt))
    meta_grp.create_dataset('Sample_taxid_ch1',data=np.array(meta_df["TaxID"], dtype=dt))
    meta_grp.create_dataset('Sample_quality',data=np.array(meta_df["QC_summary"], dtype=dt))
    meta_grp.create_dataset('Sample_title',data=np.array(meta_df["BioSample"], dtype=dt))
    meta_grp.create_dataset('Sample_type',data=np.full(len(meta_df["BioSample"]),"SRA", dtype=dt))
    h5_file.flush()
    n_srr=len(meta_df["SRR_accession"])
    n_genes=0
    genes= list()
    with bz2.BZ2File(data_file_name, 'r') as data_file:
        srr, gene_id, gene_count= data_file.readline().split(b"\t")
        first_srr = srr
        n_genes+=1
        genes.append(gene_id)
        while(True):
            srr, gene_id, gene_count= data_file.readline().split(b"\t")
            if (srr!=first_srr):
                break
            n_genes+=1
            genes.append(gene_id)
        meta_grp.create_dataset('genes',data=np.array(genes, dtype='S'))
    exp_data=data_grp.create_dataset("expression", (n_srr,n_genes),dtype= 'i4')#, compression="gzip", compression_opts=9)
    h5_file.flush()
    reader= pd.read_csv(data_file_name, sep='\t', chunksize=n_genes,names=["srr","gene","count"],iterator=True)
    srr=0
    chunk_size = 100    
    while(srr<n_srr):
        chunk = reader.get_chunk(chunk_size*n_genes)
        last_ind = (srr+chunk_size) if (srr+chunk_size)<n_srr else n_srr
        exp_data[srr: last_ind,0:n_genes]= chunk["count"].values.reshape(last_ind-srr,n_genes)
        h5_file.flush()      
        srr+=chunk_size
        if(srr%200==0):
            print(srr/n_srr)
print(1)