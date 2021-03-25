#! /usr/bin/env python

# import pandas as pd
# import sys
import datetime
# import sqlite3
from _utils import *

def download_file(dbpath="db/"):
    """
    automatically download gwas efo trait mapping files
    """
    import shutil
    import urllib.request as request
    from contextlib import closing

    url = "ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-efo-trait-mappings.tsv"
    outfile = dbpath+"/gwas-efo-trait-mappings_download_" + str(datetime.date.today()) + ".txt"

    with closing(request.urlopen(url)) as r:
        with open(outfile, 'wb') as f:
            shutil.copyfileobj(r, f)
    return(outfile)


def efo2dict():

    #if not os.path.exists(os.path.join("./db","gwas-efo-trait-mappings_download_" + str(datetime.date.today()) + ".txt")):
    #    efofile = download_file()
    #else:
    #    efofile = os.path.join("./db","gwas-efo-trait-mappings_download_" + str(datetime.date.today()) + ".txt")
    efofile = os.path.join("./db", "gwas-efo-trait-mappings.txt")
    #efofile = os.path.join("./db", "gwas-efo-trait-mappings.tsv")
    
    traits = pd.read_csv(efofile, sep='\t', names=["raw_trait","mapped_trait","efo_mappedtrait", "parent_group","efo_parent"],header=0, na_values = "NA")
    traits_sub = traits[["mapped_trait","efo_mappedtrait", "parent_group","efo_parent"]]
    traits_sub = traits_sub.drop_duplicates()

    efodict = {}
    for idx, irow in traits_sub.iterrows():
        # efodict[irow["efo_mappedtrait"].rstrip()] = [irow["mapped_trait"], irow["parent_group"], irow["efo_parent"]]
        efodict[irow["efo_mappedtrait"].rstrip()] = [irow["parent_group"], irow["efo_parent"]]

    return efodict

def add_oto(idat=None, annocol="Trait2"):

    rout = []
    efo_parent = efo2dict()

    gwasdb = query_gwasumdb()
    allefo = gwasdb.get_allefo()
    allefo_dict = dict(zip(allefo['trait'], allefo['efo']))
    # idat = idat.merge(allefo, left_on= annocol, right_on="trait")

    for idx, irow in idat.iterrows():
        itrait = irow[annocol]
        itrait_efo = "NA"

        ### first add EFO
        if itrait in allefo_dict:
            itrait_efo = allefo_dict[itrait]

        ## with EFO, add disease group
        ## using EFO, it is more accurate than mapped_trait
        if not itrait_efo or itrait_efo in ["NA", "None"]:
            rout.append(["NA", "NA", "NA"])
        else:
            if itrait_efo in efo_parent:

                rout.append([itrait_efo] + efo_parent[itrait_efo])
            else:
                if "," in itrait_efo:
                    alltrait = itrait_efo.split(",")
                    alltrait = [x.lstrip().rstrip() for x in alltrait]
                    itraits_info = []

                    for iefo in alltrait:
                        if iefo in efo_parent:
                            itraits_info.append([itrait_efo] +efo_parent[iefo])
                        else:
                            print(f"{iefo} is not in gwascatlog lib")
                            itraits_info.append(["NA", "NA", "NA"])

                    tmp = zip(*itraits_info)
                    tmp2 = [",".join(set(x)) for x in tmp]
                # tmp = pd.DataFrame(itraits_info)
                # tmp2 = tmp.agg(lambda x:",".join(x))

                    rout.append(tmp2)
                else:
                    rout.append(["NA","NA","NA"])

    return(rout)



