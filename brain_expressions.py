import pandas as pd

phenotype_annotations = pd.read_csv("annotations/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.csv")
sample_annotations = pd.read_csv("annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.csv", index_col=0)
expressions = pd.read_table("Voom_Normalized_GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct", index_col=0).T

expressions.index = map(lambda x: "-".join(x.split(".")), expressions.index)
descriptions = expressions.loc['Description']
sample_annotations = sample_annotations.loc[sample_annotations['SMTS'] == 'Brain']
merged = pd.merge(expressions, sample_annotations, how='inner',
        left_index=True, right_index=True)
merged = merged.append(descriptions).sort_index()
merged.T.to_csv("Brain_Voom_Normalized_GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct",
	sep="\t", header=True, index=True, na_rep="NA")
