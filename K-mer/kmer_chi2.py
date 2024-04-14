


import pandas as pd 
from scipy.stats import chi2_contingency
import numpy as np 
import statsmodels.stats.multicomp
import timeit
import pdb 

def kiannaChi2(a, b, fc=0.5):
    '''
    prev function call: chi2(a, b, asum, bsum, fc=0.5)
    seems like a and b are both dictionaries for a single treatment. 
    
    ex: a_treatment = {
        aaa:count,
        bbb:count,
        ...
        }

        asum/bsum of all of the counts of kmers within a single treatment group
        fc is some constant that was always set to 0.5, so I made it 0.5 unless specified. 
    '''
    asum = sum(a.values())
    bsum = sum(b.values())

    a_motifs=[]
    b_motifs=[]
    pvals=[]
    keys=list(set(list(a.keys())+list(b.keys())))

    # For kmer in kmers,
    for key in keys:
        # If kmer was found in group A, apos= num times kmer was counted in group A
        if key in a:
            apos=a[key]
        
        # If not found in group A, apos is 0
        else:
            apos=0.0
        aneg=asum-apos
        # A neg is the total counts of group A - num times kmer was counted in group A. 
        if key in b:
            bpos=b[key]
        else:
            bpos=0.0
        bneg=bsum-bpos
        if apos !=0.0 and bpos!=0.0:
            table=[[apos,aneg],[bpos,bneg]]
            #print(table)
            stat, p, dof, expected = chi2_contingency(table)
            pvals.append([key,p])
            #add bh here
    pv=list(map(list, zip(*pvals)))[1]
    sq=statsmodels.stats.multitest.multipletests(pv, method="fdr_bh")
    i=0

    # For each rejected null hyp., 
    for reject in sq[0]:
        key=pvals[i][0]
        if key in a:
            apos=a[key]
        else:
            apos=0
        # Get pval of rejected null.
        aneg=asum-apos
        if key in b:
            bpos=b[key]
        else:
            bpos=0
        bneg=bsum-bpos

        # If reject is not None
        if reject: 
            # Perform log2  enrichment
            amore=np.log2(apos/asum)-np.log2(bpos/bsum)
            bmore=np.log2(bpos/bsum)-np.log2(apos/asum)
            #print("amore: ",amore,"   bmore: ",bmore)
            if amore > fc: 
                a_motifs.append([key,sq[1][i],amore])
            elif bmore > fc: 
                b_motifs.append([key,sq[1][i], bmore])
                #motif, adjusted p-value, fold change
        i+=1
    # a_motifs are motifs enriched in treatment a and not treatment b 
    # b_motifs are motifs enriched in treatment b and not treatment a. 
        
    return a_motifs, b_motifs

def main():
    data = '/Users/matthewchaw/Desktop/School/UW/rotation_projects/Nunn/cp34h_kmer/output/2024-01-26_two_peps_all_runs_1-2-3-4.csv'
    data = pd.read_csv(data) 
    filtered_data = data[data['group'] != 'Initial']

    category = 'group'

    grouped = filtered_data.groupby(category).apply(lambda x: x.groupby('kmer')['count'].sum()).to_dict()
    kmer_dict = dict()
    for (outer_key, inner_key), value in grouped.items():
        if outer_key not in kmer_dict:
            kmer_dict[outer_key] = {}
        kmer_dict[outer_key][inner_key] = value
    start = timeit.default_timer()
    a_enriched_motifs, b_enriched_motifs = kiannaChi2(kmer_dict['Control'], kmer_dict['WCL'])

    kmer_data = {
        'group':[],
        'k-mer':[],
        'p-value':[],
        'log2 fold change':[]
    } 

    for motif in a_enriched_motifs:
        kmer_data['group'].append('enriched in control')
        kmer_data['k-mer'].append(motif[0])
        kmer_data['p-value'].append(motif[1])
        kmer_data['log2 fold change'].append(motif[2])
    for motif in b_enriched_motifs:
        kmer_data['group'].append('enriched in wcl')
        kmer_data['k-mer'].append(motif[0])
        kmer_data['p-value'].append(motif[1])
        kmer_data['log2 fold change'].append(motif[2])
    
    kmer_enrichment = pd.DataFrame(kmer_data)
    kmer_enrichment.to_csv('../output/kmer_enrichment_log2fold.csv', index=False)
    
    
    
    stop = timeit.default_timer()
    print(stop-start, 'seconds')
if __name__ == '__main__':
    main() 
