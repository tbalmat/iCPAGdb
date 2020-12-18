"""Utility functions and classes
"""

from itertools import combinations, permutations, product
from joblib import Parallel, delayed, Memory
from tqdm import trange
import stats, os, subprocess
from collections import ChainMap
import pandas as pd
import functools, time
import sqlite3, sys
from statsmodels.stats.multitest import multipletests


# import psutil
# svmem = psutil.virtual_memory()
# print("available memory: ", svmem)


def log_info(func):
    """Decorator that prints function arguments and runtime
    """
    @functools.wraps(func)
    def wrapper(args):
        print(
            "Func {} uses the following arguments:\n".format(func.__name__)
        )
        for arg in vars(args):
            print(str(arg) + "\t" + str(getattr(args, arg)))
        t1 = time.time()
        func(args)
        t2 = time.time()
        elapsed = [round(x, 2) for x in divmod(t2 - t1, 60)]
        print("\nFunction completed in  {} m {} s\n".format(elapsed[0], elapsed[1]))
    return wrapper


def check_py_version():
    """Checks the python version, exits if < 3.4."""
    py_major, py_minor = sys.version_info[:2]

    if py_major == 2:
        sys.stderr.write("CPAG requires python 3+\n")

    if py_major < 3 or py_minor < 6:
        sys.stderr.write("CPAG requires python 3 (version 3.6 or higher)\n")
        sys.exit(1)

flatten_list = lambda lst: [y for x in lst for y in x]

def intra_trait(snpdat, ncpus = 1, ldpop="EUR", cross_traits = False):
    """
    :param snpdat: a dict. Keys are trait name, and values are a list containing SNPs
    :param ncpus: int. number of threads/CPUs
    :param cross_traits: boolean. True or False.
    :return: a pandas DataFrame. the main result
    """

    print(f"Running with {ncpus} threads ...")

    if ldpop == "EUR":
        NtotalSNP = 1634900
    elif ldpop == "AFR":
        NtotalSNP = 3091723
    elif ldpop == "ASN" or ldpop == "EAS":
        NtotalSNP = 1442763

    if not cross_traits:
        assert isinstance(snpdat, dict)
        #TotalTraitNumber = len(snpdat.keys())
        #allpairs = list(combinations(range(TotalTraitNumber), 2))
        allpairs = list(combinations(snpdat.keys(),2))
        # allsnps = flatten_list(snpdat.values())
    else:
        assert isinstance(snpdat, list)
        allpairs = list(product(snpdat[0].keys(), snpdat[1].keys()))
        snpdat = dict(ChainMap(*snpdat)) ## now snpdat is new dict
        # allsnps = flatten_list(snpdat.values())

    #lddat_sub = lddat[lddat.isin(allsnps).any(axis=1)].groupby("SNP")['TAGS'].agg(set).to_dict()
    # lddb_con = query_lddb()

    # import multiprocessing as mp
    # pool = mp.Pool(processes=int(ncpus))  ## usually only use half of CPUs
    # res = [pool.apply(intra_trait_pair, args=(snpdat, allpairs[i][0], allpairs[i][1],)) for
    #        i in trange(len(allpairs), desc='Total processing: ')]
    # resall = [pool.apply_async(intra_trait_pair, args=(snpdat, allpairs[i][0], allpairs[i][1],)) for
    #           i in trange(len(allpairs), desc='Total processing: ')]
    # res = [r.get() for r in resall]

    # , require='sharedmem' , backend="multiprocessing","threading"
    res = Parallel(n_jobs=ncpus, backend="threading")(delayed(intra_trait_pair)(
        trait1snps = set(snpdat[allpairs[i][0]]), trait2snps = set(snpdat[allpairs[i][1]]),
        trait1name = allpairs[i][0], trait2name = allpairs[i][1], ldpop = ldpop)
        for i in trange(len(allpairs), desc='Total processing: '))

    res_f = pd.DataFrame(res, columns=["Trait1", "Trait2", "N1_pcut", "N2_pcut", "Nshare_direct",
                                       "Nshare_LDpairComb", "Nshare_all",  "SNPshare_all","SNPshareByLD"])

    res_f = res_f[res_f['Nshare_all'] > 0]

    if res_f.shape[0] > 0:
        res_f["P_fisher"] = res_f.apply(lambda x: cal_fisher(x,NtotalSNP), axis=1)
        res_f["N_expected"] = res_f.apply(lambda x: cal_exp(x, NtotalSNP), axis=1)
        #res_f['Padj_BH'] = multipletests(res_f["P_fisher"].astype(float), method='fdr_bh' )[1]
        res_f['Padj_FDR'] = stats.multest(res_f["P_fisher"].astype(float).tolist()).p_adj(method="BH", n=len(allpairs) )
        #res_f['Padj_Bonferroni'] = multipletests(res_f["P_fisher"].astype(float), method='bonferroni')[1]
        res_f['Padj_Bonferroni'] = stats.multest(res_f["P_fisher"].astype(float).tolist()).p_adj( method="bonferroni", n=len(allpairs) )
        res_f['Jaccard'] = res_f['Nshare_all']/(res_f['N1_pcut'] + res_f['N2_pcut'] - res_f['Nshare_all'])
        res_f['Sorensen'] = 2 * res_f['Nshare_all'] / (res_f['N1_pcut'] + res_f['N2_pcut'] )

        print("Calculating the disease/trait similarity ...")
        print("\tNote: Please use Jaccard index. The Chao-Sorensen and Chao-Jaccard methods are underpowered here. "
              "Since each SNP present just ONCE in a trait, it makes no sense to incorporate abundance!")
        res_f["ChaoSorensen"] = res_f.apply(lambda x: cal_similarity(x), axis=1)

        res_f.sort_values(by=["P_fisher"], inplace=True)
        # gwasdb = query_gwasumdb()
        # allefo = gwasdb.get_allefo() ## add efo for all traits

        # res_fo = res_f[["Trait1", "Trait2", "N1_pcut", "N2_pcut", "Nshare_direct", "Nshare_LDpairComb", "Nshare_all",
        #                 "N_expected", "P_fisher", 'Padj_BH', 'Padj_BH2', 'Padj_Bonferroni', 'Padj_Bonferroni2',
        #                 'Jaccard', 'Sorensen',
        #                 "ChaoSorensen", "SNPshare_all", "SNPshareByLD"]]

        res_fo = res_f[["Trait1", "Trait2", "N1_pcut", "N2_pcut", "Nshare_direct","Nshare_LDpairComb", "Nshare_all",
                        "N_expected", "P_fisher", 'Padj_FDR','Padj_Bonferroni', 'Jaccard', 'Sorensen', "ChaoSorensen", "SNPshare_all"]]

        return(res_fo)
    else:
        print("No overlapping traits were found, try to relax p value cutoff!")
        return(pd.DataFrame()) ## to avoid the some complain


def cal_fisher(irow, NtotalSNP):
    n1 = irow["N1_pcut"]
    n2 = irow["N2_pcut"]
    n0 = irow["Nshare_all"]
    p1 = stats.overlap_test(NtotalSNP, n1, n2, n0).fisher_t()
    return(p1)


def cal_exp(irow, NtotalSNP):
    n1 = irow["N1_pcut"]
    n2 = irow["N2_pcut"]
    n0 = irow["Nshare_all"]
    p1 = stats.overlap_test(NtotalSNP, n1, n2, n0).exp_overlap()
    return(p1)


def cal_similarity(irow):
    '''
    input is single row of pandas dataframe from CPAG output
    '''
    # dat = pd.read_csv("./output/cpag_output_H2P2pcut0.0001_20201014.out.csv", header = 0)
    # irow = dat.iloc[0]

    n1 = irow["N1_pcut"]
    n2 = irow["N2_pcut"]
    n0_all = irow["Nshare_all"]
    snplist1 = ["rs1_" + str(i) for i in range(n1 - n0_all)]
    snplist2 = ["rs2_" + str(i) for i in range(n2 - n0_all)]
    snplistO = ["rs_" + str(i) for i in range(n0_all)]
    snplist1 = snplist1 + snplistO
    snplist2 = snplist2 + snplistO
    cs_sim = stats.similarity_trait(trait1=snplist1, trait2=snplist2, method="chao_jaccard")
    return(cs_sim)


def merge_SNPpairList(l=None):
    # combine sublist from list, if the two sublist share common SNPs
    # Credit to Howards, refer to https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements
    #l = [['a', 'b', 'c'], ['b', 'd', 'e'], ['k'], ['o', 'p'], ['e', 'f'], ['p', 'a'], ['d', 'g']]

    olist = []
    while len(l)>0:
        first, *rest = l
        first = set(first)

        lf = -1
        while len(first)>lf:
            lf = len(first)

            rest2 = []
            for r in rest:
                if len(first.intersection(set(r)))>0:
                    first |= set(r) ## union and merge into first
                else:
                    rest2.append(r) ## keep it
            rest = rest2

        olist.append(first)
        l = rest

    out_pair = [list(x) for x in olist]
    n0_len = len(olist)

    return(out_pair, n0_len)


def intra_trait_pair(trait1snps = None, trait2snps = None, trait1name = None, trait2name = None, ldpop = "EUR"):
    """
    :param trait1snps: set. a set of SNPs' rsID for trait2
    :param trait2snps: set. a set of SNPs' rsID for trait2
    :param trait1name: str. trait name for trait1
    :param trait2name: str. trait name for trait2
    :return: pandas DataFrame.
    """

    lddb_con = query_lddb(ldpop = ldpop)
    lddb_con.opendb()

    # trait_i, trait_j = ["Eczema",	"Self-reported psoriasis"]
    # h2p2_prunedat = H2P2(pvcut = 1e-5)
    # snpdat = h2p2_prunedat.snpdict()
    # gwasdb = query_gwasumdb()
    # snpdat = gwasdb.search_traits(itrait= [trait1name,trait2name], pcut = 5e-100)
    # snpdat = gwasdb.search_source(source=tmptrait, pcut=5e-100)
    # print("i: {0} vs j: {1}".format(str(trait1name), str(trait2name)))
    #trait1snps = set(snpdat[trait1name])
    #trait2snps= set(snpdat[trait2name])

    n1 = len(trait1snps)
    n2 = len(trait2snps)

    SNPshare_direct = set(trait1snps) & set(trait2snps)
    n0_direct = len(SNPshare_direct)

    n1_left = trait1snps - SNPshare_direct
    n2_left = trait2snps - SNPshare_direct

    # alleftsnps = list(n1_left.union(n2_left))

    pairlist = [];    n0_extra_allldpair_len = 0

    if len(n1_left) > 0 and len(n2_left) > 0:
        #lddat1 = query_lddb(snplist=alleftsnps).groupby("SNP")['TAGS'].agg(set).to_dict()
        lddat1 = lddb_con.query_pairlist(snplist1=n1_left, snplist2=n2_left)
        # lddat1 = lddb_con.query(snplist=alleftsnps, pair=True)

        if lddat1.shape[0] > 0:
            ldpairs = lddat1.astype(str).apply(lambda x: sorted(x), axis=1)
            if len(list(ldpairs)) > 0:
                pairlist, n0_extra_allldpair_len = merge_SNPpairList(l = list(ldpairs))
            #pairlist, n0_extra_allldpair_len = dict_search(lddat1_dict, n1_left, n2_left)

    if len(pairlist) > 0:
        n0_indirect_mergepair_set = set(["-".join(x) for x in pairlist])
        shareSNPall = SNPshare_direct.union(n0_indirect_mergepair_set)
        n0_extra_len = len(n0_indirect_mergepair_set)
        n0_extra_snp = '|'.join(str(x) for x in n0_indirect_mergepair_set) # for better readable
    else:
        shareSNPall = SNPshare_direct
        n0_extra_len = 0
        n0_extra_snp = ""

    lddb_con.close()

    n0_cor = len(shareSNPall) ## n0_corrected
    # pv = stats.overlap_test(NtotalSNP, n1, n2, n0_cor)
    n0_all = "|".join(shareSNPall)
    # p_fisher = str(format(pv.fisher_t(), "10.8e"))
    # exp_overlap = str(format(pv.exp_overlap(), "10.8e"))

    out = [trait1name, trait2name, n1, n2, n0_direct, n0_extra_len, n0_cor, n0_all, n0_extra_snp]
    return(out)


def ld_clump(idat, args):
    # pop = "EUR"; #clump_r2 = 0.4; clump_kb = 10000
    # clump_p1 = 1; clump_p2 = 1; ## here no restriction of p values were applied as a SNP is significant in one dataset, but its proxy may be significant in another
    #plink_bin = "./plink_bins/plink"
    plink_bin = "plink_bins\plink.exe"
    b_file = "./db/lddat/" + args.ldpop + "_1kg_20130502_maf01"
    tag_nsnp = 10000 ## maximum number of SNPs' proxy
    #tmpfile = "tempfile4clump"
    tmpfile = "plink_bins/tmp/tempfile4clump"

    if os.path.isfile(tmpfile):
        os.remove(tmpfile)

    if not os.path.isfile(plink_bin):
        sys.exit(" {} file does not exist.".format(str(plink_bin)))

    if not (os.path.isfile(b_file + ".bed") and os.path.isfile(b_file + ".bim") and os.path.isfile(b_file + ".fam")):
        sys.exit("1KG reference {} bim/bed/fam does not exist in db/lddat folder.".format(str(b_file)))

    sys.stdout.write("Performing LD clumping for usr GWAS uploaded data ...\n\n")

    if args.snpcol not in idat.columns:
        sys.exit(" {} is not a valid column name for SNP rsID.".format(str(args.snpcol)))
    if args.pcol not in idat.columns:
        sys.exit(" {} is not a valid column name for p value.".format(str(args.pcol)))

    tmpidat = idat[[args.snpcol, args.pcol]]
    tmpidat.columns = ["SNP", "P"]
    tmpidat = tmpidat[tmpidat.SNP.notnull()]  ## remove NA and NaN
    tmpidat = tmpidat[~tmpidat["SNP"].str.contains('\.|NA') ]

    if tmpidat.shape[0] > 0:
        tmpidat.to_csv(tmpfile, sep='\t', index=False)
    else:
        sys.exit("Please use '--SNPcol' and '--Pcol' to choose columns for SNP rsID and p value.")

    cmd1 = [plink_bin, "--bfile", b_file,
            "--clump", tmpfile,
            "--clump-p1", args.ldclump_p1,
            "--clump-p2", args.ldclump_p2,
            "--clump-r2", args.ldclump_r2,
            "--clump-kb", args.ldclump_kb,
            "--out", tmpfile]
    cmd2 = " ".join([str(x) for x in cmd1])

    # prun = subprocess.Popen(cmd2, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #
    # stdout_out, stderr_out = prun.communicate()
    # if stdout_out or stderr_out:
    #     print("LD clumpping error code: ", prun.returncode)
    prun = subprocess.check_output(cmd2, shell=True, stderr=subprocess.DEVNULL)

    if os.path.exists(tmpfile + ".clumped") and os.stat(tmpfile + ".clumped").st_size > 0:
        resclump = pd.read_csv(tmpfile + ".clumped", sep='\s+')
        findat = idat[idat[args.snpcol].isin(resclump["SNP"])]
        snpdel = idat[~idat[args.snpcol].isin(resclump["SNP"])]
        resclump['SNP'].to_csv(tmpfile + ".snplist", index=False)

        if findat.shape[0] > 0:
            print("After LD clumping, " + str(findat.shape[0]) + " leading SNPs were kept")
        else:
            print("Skip LD clumping because all SNPs (rsID) are absent in 1000 genome EUR population")

        # if snpdel.shape[0] > 0:
        #     print("Removing " + str(snpdel.shape[0]) + " of " + str(idat.shape[0]) +
        #           " variants due to LD or absence from LD reference panel")
    else:
        print("Skip LD clumping because all SNPs (rsID) are absent in 1000 genome EUR population")
        findat = idat

    ldproxy_filter = {}
    # ### get LD proxy SNPs based on clumped snplist
    # if os.path.exists(tmpfile + ".clumped") and os.access(tmpfile + ".clumped", os.R_OK):
    #     cmd0 = [plink_bin, "--bfile", b_file,
    #             "--ld-snp-list", tmpfile + ".snplist",
    #             "--keep-allele-order", "--r2 in-phase with-freqs gz",
    #             "--ld-window-kb", args.ldclump_kb,
    #             "--ld-window", tag_nsnp,
    #             "--out", tmpfile + "_proxySNPs"]
    #
    #     cmd01 = " ".join([str(x) for x in cmd0])
    #     runcmd = subprocess.Popen(cmd01, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #
    #     if os.path.exists(tmpfile + "_proxySNPs.ld.gz") and os.access(tmpfile + "_proxySNPs.ld.gz", os.R_OK):
    #         ldproxy = pd.read_csv(tmpfile + "_proxySNPs.ld.gz", sep='\s+')
    #         ldproxy_cut = ldproxy[ldproxy["R2"] > args.ldclump_r2]
    #         if ldproxy_cut.shape[0] >= 1:
    #             ldproxy_filter["leadingSNPs"] = ldproxy_cut["SNP_A"].drop_duplicates().to_list()
    #             ldproxy_filter["proxySNPs"] = ldproxy_cut["SNP_B"].drop_duplicates().to_list()
    #
    #         ## in fact, just query gwas database using those leading SNPs and proxy SNPs, may reduce lots of computations
    #             total_SNPs =len(ldproxy_cut["SNP_A"].drop_duplicates()) + len(ldproxy_cut["SNP_B"].drop_duplicates())
    #             print("{} SNPs and proxys were retained in the uploaded dataset".format(total_SNPs))

    # delete temp files
    for file in os.listdir(os.getcwd()):
        if os.path.isfile(file) and file.startswith(os.path.basename(tmpfile)):
            try:
                os.remove(file)
                # continue
            except:
                print("Error while deleting file : ", file)

    return([findat, ldproxy_filter])


def check_file(dir, file):
    assert type(file) == str
    if os.path.isfile(os.path.join(dir,file)):
        return 1
    else:
        return 0
        sys.exit(" {} file is not existed.".format(str(file)))


class query_lddb(object):
    """
    create a class for h2p2 prune data
    """
    def __init__(self, ldpop= "EUR", r2cut = None):
        self.dbpath = "./db/"
        self.lddb = "cpag_gwasumstat_v1.1_" + ldpop +".ld0.4.db"
        # self.lddb = "cpag1_gwasumstat20130904." + ldpop +"_ld0.6.db"
        self.check_dbfile()
        # self._conn = self.opendb()
        self._conn = None

    def check_dbfile(self):
        if check_file(self.dbpath, self.lddb):
            return(1)
        else:
            print("{} LD sqlite3 database is missing in folder db".format(self.lddb))
            sys.exit(1)

    def opendb(self):
        conn = sqlite3.connect(self.dbpath + self.lddb, check_same_thread = False)
        conn.execute("PRAGMA synchronous = OFF")
        # conn.execute("PRAGMA journal_mode = MEMORY")
        self._conn = conn

    def close(self):
        self._conn.close()

    def query(self, snplist=None, pair=True):
        """
        :param: snplist is either single SNP or list of SNPs' rsID
        :param: pair is true, then return ldpairs within the snplist; Otherwise, return each SNPs's all LDproxy
        :return: panda dataframe
        """

        if type(snplist) == list:
            snptuple = tuple(snplist)
        elif type(snplist) == str:
            snptuple = tuple([snplist])
        else:
            print("{} is not a python list".format(snplist))
            snptuple = tuple()

        if len(snptuple) == 1:
            qry = "SELECT * FROM LDSNP WHERE SNP = '{}' OR TAGS = '{}'".format(snptuple[0], snptuple[0])
            query_out = pd.read_sql_query(qry, self._conn)
        elif len(snptuple) > 1:
            if pair:
                # return LD paris
                qry = f'''select SNP,TAGS from LDSNP where (TAGS IN {snptuple} AND SNP IN {snptuple}) AND SNP != TAGS '''
                query_out = pd.read_sql_query(qry, self._conn)
            else:
                # return SNPs' all LDproxy
                qry = f'''SELECT SNP,TAGS FROM LDSNP WHERE (SNP IN {snptuple} OR TAGS IN {snptuple}) AND SNP != TAGS  '''
                query_out = pd.read_sql_query(qry, self._conn)
        else:
            query_out = pd.DataFrame()
        #
        # if query_out.shape[0] < 1:
        #     print("No LD proxy were found in LDdb")

        return(query_out)

    def query_pairlist(self, snplist1=None,snplist2=None):
        # both snplist1 and snplist2 are set
        snptuple1 = tuple(snplist1)
        snptuple2 = tuple(snplist2)

        if len(snptuple1) == 1:
            snptuple1 = snptuple1 + ("NONE",) ## This to avoid single tuple problem
        if len(snptuple2) == 1:
            snptuple2 = snptuple2 + ("NONE",)

        qry = f''' select SNP,TAGS from LDSNP where ( (TAGS IN {snptuple1} AND SNP IN {snptuple2}) OR (TAGS IN {snptuple2} AND SNP IN {snptuple1}) ) AND SNP != TAGS    '''
        query_out = pd.read_sql_query(qry, self._conn)

        return(query_out)


class query_gwasumdb(object):
    """
    create a class for GWAS database
    """
    def __init__(self):
        self.dbpath = "./db/"
        self.gwasdb = "cpag_gwasumstat_v1.1.db"
        # self.gwasdb = "cpag1_gwasumstat20130904.db"

        self.check_dbfile()
        self._conn = self.conndb()
        self._allsource = self.get_allsource()
        self._alltraits = self.get_alltrait()

    def check_dbfile(self):
        if check_file(self.dbpath, self.gwasdb):
            return(1)
        else:
            print("{} LD sqlite3 database is not missing in folder db".format(self.lddb))
            sys.exit(1)

    def conndb(self):
        conn = sqlite3.connect(self.dbpath + self.gwasdb)
        #conn.execute("PRAGMA synchronous = OFF")
        # conn.execute("PRAGMA journal_mode = MEMORY")
        # self._conn = conn
        return(conn)

    def close(self):
        self._conn.close()

    def get_allsource(self):
        qry = "SELECT distinct(source) FROM GWAStable"
        allsource = pd.read_sql_query(qry, self._conn)

        return(list(allsource["source"]))

    def get_alltrait(self):
        qry = "SELECT distinct(trait) FROM GWAStable"
        alltraits = pd.read_sql_query(qry, self._conn)

        return(list(alltraits["trait"]))

    def get_allefo(self):
        qry = "SELECT DISTINCT trait, efo FROM GWAStable"
        allefo = pd.read_sql_query(qry, self._conn)
        return (allefo)

    def search_source(self, source = None, pcut= 5e-8):

        if source in ["All","all"]:
            source = self.get_allsource()

        if type(source) == str:
            source = [source]
        assert isinstance(source, list)

        sourcetuple = tuple(source)

        if len(set(source) & set(self._allsource)) < 1:
            print("{} is not a valid GWAS data source, please choose from {}".format(source, self._allsource))

        if len(source) == 1:
            qry = "SELECT distinct trait,SNP FROM GWAStable WHERE source = '{}' AND pval <= '{}' ".format(sourcetuple[0], pcut)
            res = pd.read_sql_query(qry, self._conn)
        elif len(source) > 1:
            qry = f"SELECT distinct trait,SNP FROM GWAStable WHERE source IN {sourcetuple} AND pval <= {pcut} "
            res = pd.read_sql_query(qry, self._conn)

        #res = res.drop_duplicates()
        if res.shape[0] < 1:
            print("No SNPs were left with pvalue < {0} for {1}".format(pcut, source))
        else:
            print("{0} traits and {1} unique SNPs were kept for pvalue <= {2} for {3}".format(len(set(res["trait"])), len(set(res["SNP"])), pcut, source))

        resdict = {k: list(set(v)) for k, v in res.groupby("trait")["SNP"]}

        return(resdict)

    def search_traits(self, itrait = None, pcut= 5e-8):

        if type(itrait) == str:
            if len(set([itrait]) & set(self._alltraits)) < 1:
                print(" '{}' traits/diseases are not found in CPAG db v1.0!".format(itrait))

            qry = "SELECT distinct trait,SNP FROM GWAStable WHERE trait = '{}' AND pval < '{}' ".format(itrait, pcut)
            res = pd.read_sql_query(qry, self._conn)
        else:
            if type(itrait) == list:
                if len(set(itrait) & set(self._alltraits)) < 1:
                    print(" '{}' traits/diseases are not found in CPAG db v1.0!".format(itrait))

                qry = f"SELECT distinct trait,SNP FROM GWAStable WHERE trait IN {tuple(itrait)} AND pval < {pcut} "
                res = pd.read_sql_query(qry, self._conn)
            else:
                print(" '{}' is not a valid GWAS traits/disease!".format(itrait))
                res = pd.DataFrame()

        resdicts = {k: list(set(v)) for k, v in res.groupby("trait")["SNP"]}

        return(resdicts)
