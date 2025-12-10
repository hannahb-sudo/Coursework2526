#!/usr/bin/ python 3
"""
Title: Python script for Gene Expression Analysis
Author: Hannah Byrne
Date: 11/12/25
"""
#Import Libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


#Read in TSV files of each gene comparison
Dataset1= pd.read_csv("C:/Users/hanna/OneDrive/Documents/Computational Bio/Coding Challenge/Python  coursework/Gene_Expression/A_vs_B.deseq2_results.tsv", sep ="\t")
Dataset2 = pd.read_csv("C:/Users/hanna/OneDrive/Documents/Computational Bio/Coding Challenge/Python  coursework/Gene_Expression/A_vs_E.deseq2_results.tsv", sep ="\t")


"""
Data Quality Checks - 
Missing Values and Duplicates
"""
#Function to check for missing values
def missing_values(deseq):
    missing_values = deseq.isnull().sum()
    print("Missing values: ", missing_values)

#Function to check for duplicates
def duplicate_rows(deseq):
    duplicate_rows = deseq[deseq.duplicated()]
    print("Duplicate rows: ", duplicate_rows)


"""
Summary Statistics -
The Number of Significantly Up and Down Regulated Genes per Dataset
Statistics Table for the P-values and log2FC per Dataset
"""
#Function for finding the significantly up and down regulated genes
def signif_gene_number(deseq):
    d = deseq.dropna() #remove NA
    upreg = d.query("pvalue < 0.05 & log2FoldChange >= 1") #filter for genes with significant P value and significant pvalue
    upreg_total = len(upreg.index) #total number of genes that fit parameters
    print(upreg_total, "upregulated genes")
    downreg = d.query("pvalue < 0.05 & log2FoldChange <= -1 ") #filter for genes with negative logFC and significant pvalue
    downreg_total = len(downreg.index) #total number of genes that fit parameters
    print(downreg_total, "downregulated genes")

#Filter for just p-values and logFC and do summary stat table
def summary(deseq):
    d = deseq.dropna()
    filtered_deseq = d.filter(["pvalue", "log2FoldChange"])
    print(filtered_deseq.describe())

"""
Plots of Gene Expression -
Volcano Plot Showing significance and Magnitude of Changes per Dataset
MA Plot displaying the Relationship between log2FC and Mean Expression per Dataset
Histogram of P-values Showing Distribution of Statistical Significance per Dataset
Heatmap of the Top Differentially Expressed Genes across the Datasets
"""
#Volcano Plot Function
def volcano_plot(deseq):
    plt.scatter(x=deseq["log2FoldChange"], y=deseq["padj"].apply(lambda x: -np.log10(x)), alpha=0.5, s=2, color="grey",
                label="Not Significant") #Plot all genes

    down = deseq[(deseq["log2FoldChange"] <= -1) & (deseq["padj"] <= 0.01)] #Define downregulated genes
    up = deseq[(deseq["log2FoldChange"] >= 1) & (deseq["padj"] <= 0.01)] #Define upregulated genes

    plt.scatter(x=down["log2FoldChange"], y=down["padj"].apply(lambda x: -np.log10(x)), s=4, color="blue",
                label="Downregulated") #Plot downregulated genes
    plt.scatter(x=up["log2FoldChange"], y=up["padj"].apply(lambda x: -np.log10(x)), s=4, color="hotpink",
                label="Upregulated") #Plot upregulated genes

    plt.xlabel("Log2FoldChange")
    plt.ylabel("-log10(p-adj)")
    plt.axvline(x=-1, color="grey", linestyle="--") #Lines corresponding to the expression thresholds set
    plt.axvline(x=1, color="grey", linestyle="--")
    plt.axhline(y=(0.05), color="grey", linestyle="--")
    plt.legend()
    plt.title("Volcano Plot of Significant Upregulated and Downregulated Genes")
    plt.show()

#MA Plot Function
def MA_plot(deseq):
    plt.xscale("log")
    plt.scatter(x=deseq["baseMean"], y=deseq["log2FoldChange"], alpha=0.5, s=1, color="lightseagreen",
                label="Not Significant") #Plot all genes

    down = deseq[(deseq["baseMean"] <= 10000) & (deseq["log2FoldChange"] < -1)]  # Define downregulated genes
    up = deseq[(deseq["baseMean"] <= 10000) & (deseq["log2FoldChange"] > 1)] #Define upregulated genes

    plt.scatter(x=down["baseMean"], y=down["log2FoldChange"], s=1, color="darkorange",
                label="Downregulated")  # Plot downregulated genes
    plt.scatter(x=up["baseMean"], y=up["log2FoldChange"], s=1, color="hotpink", label="Upregulated") #Plot upregulated genes
    plt.xlabel("log(Base Mean)")
    plt.ylabel("log2FC")
    plt.axhline(y=0, color="black", linewidth=0.5) #Add in y = 0 line
    plt.legend()
    plt.title("MA Plot of logFC and BaseMean")
    plt.show()

#Histogram Function
def histo_gene(deseq):
    plt.hist(deseq["pvalue"], bins=50, color="hotpink")
    plt.xlabel("P-value")
    plt.ylabel("Frequency")
    plt.title("Histogram of P-value")
    plt.show()

#Heatmap Function
def merge_heatmap(data1, data2):
    data1["dataset"] = "B" #Add column identifying dataset origin
    data2["dataset"] = "E"
    comps = pd.concat([data1, data2]) #concatenate new dataframes
    comps_filtered = comps.dropna()  # Remove missing values
    signif_comps = comps_filtered.query("pvalue < 0.05")  # Significant differentially expressed
    # Identify top differentially expressed genes
    top_upreg = signif_comps.query("log2FoldChange > 15")
    top_downreg = signif_comps.query("log2FoldChange < -15")
    top_diff = pd.concat([top_upreg, top_downreg])
    # Select all top differentially expressed gene IDs present in concatenated dataset
    diff_gene = comps[comps["gene_id"].isin(top_diff["gene_id"])]
    heat = diff_gene.pivot(index="gene_id", columns="dataset", values="log2FoldChange") #pivot dataframe
    sns.heatmap(heat, cmap="crest").set(title="Heatmap of Top Differentially Expressed Genes") #Seaborn heatmap
    plt.show()

"""
Significant Genes List - 
Significanty Upregulated and Downregulated Genes with Their Corrsponding log2FC, P-value and Adjusted P-values
"""
#Significant Gene List Function
def signif_gene_list(data1, data2):
    data1["dataset"] = "B"  # Add column identifying dataset origin
    data2["dataset"] = "E"
    comps = pd.concat([data1, data2])  # concatenate new dataframes
    comps_filtered = comps.dropna()  # Remove missing values
    signif_comps = comps_filtered.query("pvalue < 0.05") # Filter for significant pvalue
    signif_upreg = signif_comps.query("log2FoldChange >= 1") #Positive log2FC value for upregulated genes
    signif_downreg = signif_comps.query("log2FoldChange <= -1") # Negative log2FC value for downregulated genes
    signif_reg = pd.concat([signif_upreg, signif_downreg])
    signif_reg.filter(items=["gene_id", "log2FoldChange", "pavlue", "padj", "dataset"])
    print(signif_reg)


"""
Run Functions on datasets
"""
def main():
    print("Data checks")
    print("Missing Values:")
    print("Dataset 1")
    missing_values(Dataset1)
    print("Dataset 2")
    missing_values(Dataset2)

    print("Duplicate Rows:")
    print("Dataset 1")
    duplicate_rows(Dataset1)
    print("Dataset 2")
    duplicate_rows(Dataset2)

    print("Summary Stats:")
    print("The Number of Significantly Up and Down Regulated Genes per Dataset")
    print("Dataset 1")
    signif_gene_number(Dataset1)
    print("Dataset 2")
    signif_gene_number(Dataset2)
    print("Statistics Table for the P-values and log2FC per Dataset")
    print("Dataset 1")
    summary(Dataset1)
    print("Dataset 2")
    summary(Dataset2)

    volcano_plot(Dataset1)
    volcano_plot(Dataset2)
    MA_plot(Dataset1)
    MA_plot(Dataset2)
    histo_gene(Dataset1)
    histo_gene(Dataset2)
    merge_heatmap(Dataset1, Dataset2)
    print("Significant Genes List")
    signif_gene_list(Dataset1, Dataset2)
    pass

if __name__ == "__main__":
    main()


