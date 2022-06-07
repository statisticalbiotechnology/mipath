import pandas as pd
import numpy as np
import os.path
import tarfile
import requests

cache_path = "data/reactome"

def get_cache_path():
    """ Return a path suitable to store cached files """
    try:
        os.mkdir(cache_path)
    except FileExistsError:
        pass
    return cache_path

class FileCache:
    """
    Object handling file names, mostly an abstract base class
    """
    def __init__(self, file_name):
        self.file_name = file_name

    def get_file_name(self):
        return self.file_name

class UrlFileCache(FileCache):
    """
    Object handling files with urls wich are downloaded and cashed locally
    """
    def __init__(self, file_name, url):
        self.file_name = file_name
        self.url = url

    def track_dl(self):
        response = requests.get(self.url, stream=True)
        with open(self.file_name, "wb") as handle:
            for data in response.iter_content():
                handle.write(data)

    def get_file_name(self):
        if not os.path.isfile(self.file_name):
            self.track_dl()
        return self.file_name

class TsvFileTracker(FileCache):
    """
    Object handling tab separated files
    """
    def __init__(self, file_name, filler):
        self.file_name = file_name
        self.filler = filler

    def get_file_name(self):
        return self.file_name

    def read_file(self):
        if not os.path.isfile(self.get_file_name()):
            df = self.filler()
            self.write_file(df)
            return df
        return pd.read_csv(self.get_file_name(), sep="\t", index_col=0)

    def write_file(self,df):
        df.to_csv(self.file_name, sep="\t")

def download_file(path, url):
    "This function downloads a file, path, from an url, if the file is not already cached"
    if not os.path.isfile(path):
        response = requests.get(url, stream=True)
        with open(path, "wb") as handle:
            for data in response.iter_content():
                handle.write(data)
    return path

def get_reactome_raw(organism = "HSA", gene_anot = "Ensembl"):
    reactome_url = "https://reactome.org/download/current/"
    reactome_fn = "2Reactome_All_Levels.txt"
    fn = gene_anot + reactome_fn
    path = os.path.join(get_cache_path(),fn)
    url = reactome_url + fn
    reactome_df = pd.read_csv(download_file(path, url),
                        sep='\t',
                        header=None,
                        usecols=[0,1,3],
                        names=["gene","reactome_id","reactome_name"])
    organism = "R-" + organism
    reactome_df = reactome_df[reactome_df["reactome_id"].str.startswith(organism) ]
    return reactome_df


def get_reactome(organism = 'HSA', gene_anot = 'Ensembl'):
    reactome_raw = get_reactome_raw(organism, gene_anot)
    reactome = pd.concat((reactome_raw.groupby('reactome_id')['gene'].unique(), reactome_raw.groupby('reactome_id')['reactome_name'].first()), axis=1)
    reactome.columns = ['genes','name']
    return reactome
    
    
def get_pathway(data, reactome, pathway, min_genes = 1):
    genes_path = reactome.loc[pathway, 'genes']
    genes_intesect = np.intersect1d(genes_path, data.columns)
    if len(genes_intesect) < min_genes:
        raise ValueError('Not enough genes in pathway/data')
    out_df = data[genes_intesect]
    return out_df

def get_reactome_hugo(hugo_dic, organism="HSA",):
    reactome_raw = get_reactome_raw(organism=organism, gene_anot="Ensembl")
    hugo_df = pd.read_csv(hugo_dic, sep = '\t', index_col=1)
    reactome_raw = reactome_raw[reactome_raw['gene'].isin(hugo_df.index)]
    reactome_raw['gene'] = hugo_df.loc[reactome_raw['gene']].values
    reactome = pd.concat((reactome_raw.groupby('reactome_id')['gene'].unique(), reactome_raw.groupby('reactome_id')['reactome_name'].first()), axis=1)
    reactome.columns = ['genes','name']
    return reactome

def parse_gmt(gmt_path, gene_id_dict_path = None):
    #gene_id_dict_path is the path to a file exported from biomart with only Gene Stable ID and Gene Name fields
    if gene_id_dict_path is not None:
        dict = pd.read_csv(gene_id_dict_path, sep='\t', index_col=1).dropna()
    gmt = pd.DataFrame(columns=['genes', 'url'])
    with open(gmt_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.split('\t')
            gmt.loc[line[0],'url'] = line[1]
            if gene_id_dict_path is not None:
                genes = dict[dict.index.isin(line[2:])].values.T.tolist()[0]
            else:
                genes = line[2:]
            gmt.loc[line[0],'genes'] = genes
    return gmt
