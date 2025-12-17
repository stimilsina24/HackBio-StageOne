##--------------------------------------------------------------------------------------------------------------------------------------------------------------##
##Function to perform DNA to Protein translation, based on user inputed DNA##

def translate(): #define the function for protein translation

    #Dictionary matching codon sequence with amino acid
    codon_table = {
'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
'TAT': 'Y', 'TAC': 'Y', 'TAA': 'STOP', 'TAG': 'STOP', # Stop codon
'TGT': 'C', 'TGC': 'C', 'TGA': 'STOP', 'TGG': 'W', # Stop codon
'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', # Start codon
'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}
    DNA = input('DNA sequence:')
    protein_sequence = "" #empty protein sequence variable to add the translated protein
    start_codon = DNA.find("ATG") #Define start codon

    #If statement to skip empty or length <3 sequences, t erminate function
    if len(DNA) < 3: #skip sequences with less than 3 bases
        print("Warning: sequence length less than 3..terminating")
        return #terminate the function if length less than 3

    #if statement to skip invalid sequence string
    ok_bases = {'A','T','G','C'} #creates a set consisting of bases allowed in sequence
    if not set(DNA).issubset(ok_bases): #divides DNA into a set of individual characters then compares for overlapping characters with ok_bases set
        invalid = set(DNA) - ok_bases #subtract 2 sets of characters to get invalid bases/characters that are not in ok_bases
        print(f"Warning: Non-DNA character {invalid} found in sequence..skipping")
        return #terminate the function

    #for statement to run the function once 
    for i in range(start_codon, len(DNA)-2,3): 
        codon = DNA[i:i+3]
        if codon in codon_table:
            amino_acid = codon_table[codon]
            if amino_acid == 'STOP':
                break
            protein_sequence += amino_acid

    print(f"Protein sequence: {protein_sequence}")
    
translate() #Run the function

#Sample DNA sequence for translate input
#slack = "Santosh T"
#X = "Datanewb2"

##--------------------------------------------------------------------------------------------------------------------------------------------------------------##

##Hamming Distance Function for Slack and X usernames, based on user inputed strings##
#Import packages
from scipy.spatial import distance

##Function to perform hamming distance
def hamming_dist():
    char_slack = input('Slack username:')
    char_X = input('X username:')
    if len(char_slack) != len(char_X):
        print("Make sure the lengths of the two strings are the same")
        return
    hamming_dist = 0
    for i in range(len(char_slack)):
        if char_slack[i] != char_X[i]:
            hamming_dist += 1
    print(f"Hamming distance between your slack and X handle is {hamming_dist}")
    
hamming_dist()

#Sample strings for input into the hamming distance function
#slack = "Santosh T"
#X = "Datanewb2"

##--------------------------------------------------------------------------------------------------------------------------------------------------------------##
##Figure panel assignment##

#Import all packages
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.image as mpimg



#Set current working directory and figures directory
figures = os.path.join(os.getcwd(), 'figures') #Define figures directory
figures
os.getcwd()

##Part A- Gene Expression Analysis (a-b)##

##Import gene expression datasets-

##Normalized counts for HBR vs UHR samples
counts_url = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv"

## Differential expression results (chromosome 22)
de_url = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv"

#Import csv as dataframes
df_counts = pd.read_csv(counts_url, index_col = 0) #Normalized counts
df_de = pd.read_csv(de_url, index_col = 0) # Chromosome 22, DE results

df_counts.head() #Preview the counts dataset
df_de.head() #Preview the differential expression dataset

#Clustermap of data h
cmap_fig = 'Clustermap_fig.png'
cmap_fig_path = os.path.join(figures, cmap_fig) # Join figure directory and filename
cmap = sns.clustermap(df_counts, cmap = "Blues",linewidths = 1,linecolor = "Black", figsize=(4,8))
cmap.savefig(cmap_fig_path, dpi=300, bbox_inches='tight')

#Create a figure with subplots using matplotlib
fig, axes = plt.subplots(2,3,figsize = (12,8)) #create subplots with 2 rows, 3 columns and size of 12 x 8
axes = axes.flatten() #Flatten axes to 1d array for easier access

## a) Heatmap generation - clustered heatmap of the top differentially expressed genes between HBR and UHR samples.
##Label both genes and samples. 
##Use a color gradient (e.g., Blues) to indicate expression levels.

#Define directory for heatmap figure
hm_fn = 'Figure-panel.png' #define filename
hm_path_fig = os.path.join(figures, hm_fn) # Join figure directory and filename

hm = sns.clustermap(df_counts, cmap = "Blues",linewidths = 1,linecolor = "Black", figsize=(4,8))
hm.fig.canvas.draw()
hm.savefig(hm_path_fig, dpi=300, bbox_inches="tight") #Save the figure in the 
plt.close(hm.fig)

# plt.close(hm.fig)
img = plt.imread(hm_path_fig)
axes[0].imshow(img, aspect = 0.65) #add the heatmap on 1st subplot, extend to fill the subplot space
axes[0].set_title("a", loc = "left", x = -0.75) #label subplot as a
axes[0].axis("off")
plt.tight_layout()

## b) Volcano plot
df_de.head() #Preview the differential expression dataset

##Dictionary to define a custom color palette for each group
volcano_colors = {
    "ns": "grey", #color for not significant genes
    "up": "green", #color for upregulated genes
    "down": "orange" #color for downregulated genes
}

v = sns.scatterplot(df_de, x = 'log2FoldChange', y = '-log10PAdj', hue = "significance", palette = volcano_colors, ax = axes[1]) #scatterplot with custom color palette
v.axhline(y = 0, linestyle = '--') #horizontal dashed line at y =0
v.axvline(x = -1, linestyle = '--') #vertical dashed line at x =-1
v.axvline(x = 1, linestyle = '--') #horizontal dashed line at x =1
v = axes[1].set_title("b", loc = "left", x = -0.2)

##Part B- Breast Cancer Data Exploration (c-f)##

##Import breast cancer diagnostic data-

##Breast Cancer Wisconsin dataset
brca_url = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"

df_brca = pd.read_csv(brca_url, index_col = 0)
df_brca.head()
print(df_brca.info())

## c) Scatter plot(radius vs texture
sns.scatterplot(data = df_brca, x = 'radius_mean', y = 'texture_mean', hue = 'diagnosis', ax = axes[2])
axes[2].set_title("c", loc = "left", x = -0.2)

## d) Correlation Heatmap
df_brca_subset = df_brca.loc[:,['radius_mean', 'texture_mean', 'perimeter_mean', 'area_mean', 'smoothness_mean', 'compactness_mean']]
correlation_mat = df_brca_subset.corr()

sns.heatmap(correlation_mat, annot = True, ax = axes[3], cmap = "Blues")
axes[3].set_title("d", loc = "left", x = -0.5)

# print(df_brca_subset)

## e) Scatter Plot(smoothness vs compactness)
sns.scatterplot(data = df_brca, x = 'smoothness_mean', y = 'compactness_mean', hue = 'diagnosis', ax = axes[4])
axes[4].set_title("e", loc = "left", x = -0.2)

## f) Density plot(area distribution)
sns.kdeplot(data = df_brca, x = 'area_mean', hue = 'diagnosis', fill = 'True', ax = axes[5])
axes[5].set_title("f", loc = "left", x = -0.2)

fig.tight_layout() #Prevent overlap in the figure with panels

#Save the figure
panel = 'Figure-panel.png' #define filename
panel_path = os.path.join(figures, panel) # Join figure directory and filename
plt.savefig(panel_path)
plt.show() #Print the figure
# plt.close()

