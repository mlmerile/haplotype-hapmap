import numpy as np
import pandas as pd

def create_df_family_links(filename, l_pop_selected=None):
    cols = ['FID', 'IID', 'dad', 'mom', 'sex', 'pheno', 'population']
    df_family_link = pd.read_csv(filename, sep='\t', names=cols)

    # Select the populations
    if l_pop_selected != None:
        df_family_link = df_family_link.loc[
            df_family_link['population'].isin(l_pop_selected), :
        ]

    return df_family_link

def filter_trio(df):
    condition = ((df['IID'] != '0') &
                 (df['dad'] != '0') &
                 (df['mom'] != '0'))
    return df.loc[condition,:]

def hapmap_to_df(filenames, family_links):

    l_df_phased = []

    for f in filenames:
        df_haplo_phased_pop = pd.read_csv(f, sep='[\t,\s+]')
        df_haplo_phased_pop = df_haplo_phased_pop.transpose().ix[2:-1]
        l_df_phased.append(df_haplo_phased_pop)

    df_haplo_phased = pd.concat(l_df_phased)

    df_haplo_childs = pd.DataFrame(columns=df_haplo_phased.columns,
                                   index=range(len(df_haplo_phased.index)/2))

    j=0
    for i, ind in family_links.iterrows():
        dad = ind['dad']
        mom = ind['mom']
        pop = ind['population']

        if dad != '0' and mom != '0':
            dad_label = dad + '_A'
            mom_label = mom + '_A'
            if (dad_label in df_haplo_phased.index and
                mom_label in df_haplo_phased.index):

                df_haplo_childs.loc[j] = (df_haplo_phased.loc[dad_label])
                df_haplo_childs.loc[j+1] = (df_haplo_phased.loc[mom_label])
                j += 2

            else:
                print("Dad or Mom not found in the dataset")
        else:
            print("Not a Trio")
                

    df_haplo_childs.dropna()

    return df_haplo_childs
                
def nucleo_to_01(df):
    def convert(series):
        series_most_frequent = series.value_counts().idxmax()
        series.loc[(series == series_most_frequent)] = 0
        series.loc[(series != 0)] = 1
        return series

    df.apply(convert)
    df = df.astype(np.int8)
    return df

def build_array(filenames, family_links_filenames, l_pop):

    fam = filter_trio(create_df_family_links(family_links_filenames, l_pop))

    df_haplotypes = nucleo_to_01(hapmap_to_df(filenames, fam))

    mat_haplotypes = df_haplotypes.as_matrix()
    mat_genotypes = mat_haplotypes[::2,:] + mat_haplotypes[1::2,:]

    return (mat_genotypes, mat_haplotypes)