import pandas as pd
data = pd.read_csv('/cluster/projects/novartis_antibody/csv_data/may_2017/HiSeq/CDR3_counts.csv')
data_enbrel_a = data[['Unnamed: 0', 'Lucentis_a_R3']][(data['Lucentis_a_R2']>0) | (data['Lucentis_a_R3']>0)]
#data_enbrel_a = data[['Unnamed: 0', 'Lucentis_a_R3']][(data['Lucentis_a_R1']>0) | (data['Lucentis_a_R2']>0) | (data['Lucentis_a_R3']>0)]
data_enbrel_a.to_csv('dataset/lucentis_a_R2_R3.tsv', sep='\t', header=False, index=False)
