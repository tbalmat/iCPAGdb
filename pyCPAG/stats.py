from __future__ import division
import math, random
import scipy, sys
import scipy.stats
import numpy as np
import pandas as pd
    
class multest():
# ## The multiple testing correction was from
    def __init__(self, p):
        self._plist = p

    def p_adj(self, method, n=0):
      import numpy as np

      if not method:
        method = 'BH'
      elif method == "fdr" or method == "FDR":
        method = "BH"
      elif method not in ["fdr","FDR", "bonferroni", "hochberg", "BH"]:
        print("Please use a vaild correction method for p value , available methods are ",  ["fdr","FDR", "bonferroni", "hochberg", "BH"])

      p = np.array(self._plist)
      if not n:
          n = len(p)
          
      if n < len(p):
          print("{} must be larger than length of pvalues".format(n))
          sys.exit(1)
          
      p0 = p
      nna = np.isnan(p) ^ True

      p = p[nna]
      lp = len(p)

      if n <= 1:
        return(p0)

      if method == 'bonferroni':
        p0[nna] = np.minimum(1, n * p)

      elif method == 'hochberg':
        i = np.arange(lp, 0, -1.)
        o = np.argsort(p)
        o = o[::-1]
        ro = np.argsort(o)
        p0[nna] = np.minimum(1, np.minimum.accumulate((n - i + 1.) * p[o] ))[ro]

      elif method == 'BH':
        ## This is Benjamini-Hochberg
        i = np.arange(lp, 0, -1.)
        o = np.argsort(p)
        o = o[::-1]
        ro = np.argsort(o)
        p0[nna] = np.minimum(1, np.minimum.accumulate(n / i * p[o]))[ro]

      return(p0)


class overlap_test(object):
    # __slots__ = ('_n1','_n2','_ntotal', '_overlap', '_pzero_overlap', 'pr_overlap', 'pr_notoverlap')
    def __init__(self, ntotal, n1, n2, overlap):
        if ntotal == 0:
            ntotal = 1
        self._ntotal = ntotal
        self._n1 = n1
        self._n2 = n2
        self._overlap = overlap
        self._pzero_overlap = math.pow(1- self._n1 * self._n2/math.pow(self._ntotal,2), self._ntotal)
        self.pr_overlap = math.pow(self._n1 * self._n2/math.pow(self._ntotal,2), self._overlap)
        self.pr_notoverlap = math.pow(1- self._n1 * self._n2/math.pow(self._ntotal, 2), self._ntotal - self._overlap)

        if self._overlap > self._n1 or self._overlap > self._n2 :  # some times, the shared snps is larger than n1 or n2, because we count by genes, not by SNPs, then self._overlap equal the minimum of (self._n1, self._n2)
            self._overlap = min(self._n1, self._n2)

    def fisher_t(self):
        """Fisher's test takes the 2x2 table
                a b
                c d
            and interprets the situation under the null hypothesis as taking a
            sample of size a+c from an urn with a+b white balls and c+d black balls."""
        n = [[self._overlap, self._n1 - self._overlap], [self._n2 - self._overlap, self._ntotal - self._n1 - self._n2 + self._overlap ]]
        # n = [[4, 48 - 4], [45-4, 1634900 - 48 - 45 + 4]]
        oddratio, p = scipy.stats.fisher_exact(n, alternative="greater")
        return(p)

    def hyperm_pmf(self):
        """Pr(X=k) """
        hyp = scipy.stats.hypergeom.pmf(self._overlap, self._ntotal, self._n1, self._n2)   ## n1 markered, Sampling n2
        return(hyp)

    def hyperm_sf(self):
        """ This is hypergeometic test
         Pr(X > k) GREATER than P(X = k) """
        psf = scipy.stats.hypergeom.sf(self._overlap-1, self._ntotal, self._n1, self._n2)
        return(psf)

    def hyperm_cdf(self):
        """ Pr(X <= k)   LESS than P(X = k)"""
        psf = scipy.stats.hypergeom.cdf(self._overlap-1, self._ntotal, self._n1, self._n2)

        return(psf)

    def biomial_p(self):
        """  binomal proportion test
            Pr(X > k) = 1 - sum_Pr(X <=K)
           """
        if self._n1 == self._ntotal:
            self._ntotal = self._n1 + 1
        proportion = (self._n2 - self._overlap)/( self._ntotal - self._n1)
        psf = scipy.stats.binom.sf(self._overlap, self._n1, proportion)
##        exp = scipy.stats.binom.expect(self._overlap, self._n1, proportion)
        return(psf)

    def ks_2sample(self):
        genset1 = []
        genset2 = []

        k = scipy.stats.ks_2samp(genset1, genset2)
        return(round(k,4))

    def exp_overlap(self):
        """ From wiki, the hypergeometic expect value is n*K/N"""
        exp= self._n1 * self._n2/self._ntotal
        return(exp)

    def var_overlap(self):
        """ From wiki, the hypergeometic expect variance"""
        v = self._n1 * self._n2 * (self._ntotal - self._n1) * (self._ntotal - self._n2)/(self._ntotal/self._ntotal - 1)
        return(v)


    def simulation(self, nperm): ## Monte Carolo Simulation
        cn = 0
        for i in xrange(int(nperm)):
            list1 = set(random.sample(xrange(self._ntotal), self._n1))
            list2 = set(random.sample(xrange(self._ntotal), self._n2))
            op = list1.intersection(list2)
            if len(op) > self._overlap:
                cn += 1

        return(round(cn/float(nperm),4))


def make_matrix(trait1, trait2):
    # trait1 = ["rs12", "rs13", "rs14", "rs15"]
    # trait2 = ["rs12", "rs15", "rs16", "rs17", "rs19"]
    alltrait = list(set(trait1).union(set(trait2)))

    df = np.zeros((2,len(alltrait)), dtype=int)
    check_ele = lambda x,trait: 1 if x in trait else 0

    for i in range(len(alltrait)):
        df[0,i] = check_ele(alltrait[i], trait1)
        df[1,i] = check_ele(alltrait[i], trait2)

    outdf = pd.DataFrame(df, index = ["trait1","trait2"], columns=alltrait)
    return(outdf)


def similarity_trait(trait1=None, trait2 = None, method = "chao_sorensen",version = "rare", freq="None"):
    '''
    :param trait1
    :param trait2
    :param method, select method from jaccard, sorensen, chao-jaccard, chao-sorensen
    :param version, only for chao-sorensen or chao-jacard, choose one from [prob or rare].
    :param freq, a vector contrains frequency, only for chao-sorensen

    :return similarity value, a float
    '''

    # trait1 = ["rs12", "rs13", "rs14", "rs15"]
    # trait2 = ["rs12", "rs15", "rs16", "rs17","rs19"]

    avail_method = ["jaccard", "chao_jaccard","chao_sorensen", "sorensen", "cosine","simpson","geometric"]

    sim = 0

    if method == "jaccard":
        common_snp = set(trait2).intersection(set(trait1))
        C = len(common_snp)
        A = len(trait1)
        B = len(trait2)
        sim = C/(A+B-C)

    elif method == "sorensen":
        common_snp = set(trait2).intersection(set(trait1))
        C = len(common_snp)
        A = len(trait1)
        B = len(trait2)
        sim = 2 * C / (A + B)

    elif method == "cosine":
        common_snp = set(trait2).intersection(set(trait1))
        C = len(common_snp)
        A = len(trait1)
        B = len(trait2)
        if A >0 and B > 0:
            sim = C / np.sqrt(A * B)
        else:
            sim = 0

        if sim > 1:
            sim = 1

    elif method == "simpson":
        common_snp = set(trait2).intersection(set(trait1))
        C = len(common_snp)
        A = len(trait1)
        B = len(trait2)
        if A >0 and B > 0:
            sim = C / np.min([A,B])

        if sim > 1:
            sim = 1

    elif method == "geometric":
        common_snp = set(trait2).intersection(set(trait1))
        C = len(common_snp)
        A = len(trait1)
        B = len(trait2)
        if A >0 and B > 0:
            sim = C ** 2 /(A * B)

        if sim > 1:
            sim = 1

    elif method in ["chao_sorensen", "chao_jaccard"]:
        # print("Chao-Sorensen and Chao-Jaccard methods are modified and not abundance-based, "
        #       "since each SNP is assumed to occur just ONCE in trait1/trait2 ")

        if version == "prob" and freq == "None":
            print("The prob version of the Chao index is not available for incidence-based data."
                  "Argument freq was ignored and thus the abundance-based formulae was used.")

        dat = make_matrix(trait1, trait2)
        colfilter = dat.sum(axis = 0) > 0
        rowfilter = dat.sum(axis = 1) > 0
        datn = dat.loc[rowfilter,colfilter]

        sim0 = np.zeros((datn.shape[0], datn.shape[0]), dtype=float)

        for i in range(dat.shape[0]-1):
            samp1 = dat.iloc[i,:]
            for j in range(i+1,dat.shape[0]):
                samp2 = dat.iloc[j,]

                pair = np.stack((samp1,samp2))
                pair_pa = pair.copy()
                pair_pa[pair_pa>0] = 1 ## exist and absent
                D12_TF = pair_pa.sum(axis=0) == 2

                if np.all(D12_TF == False):
                    sim0[j,i] = 0
                else:
                    D12_spp = pair[:,D12_TF]
                    n = pair.sum(axis=1)[0]
                    m = pair.sum(axis=1)[1]
                    if version == "prob":
                        U = np.sum(D12_spp[0,]/n)
                        V = np.sum(D12_spp[1,]/m)
                    if version == "rare":
                        f1_ = np.sum(D12_spp[0,] == 1)
                        f_1 = np.sum(D12_spp[1,] == 1)
                        f2_ = np.sum(D12_spp[0,] == 2)
                        f_2 = np.sum(D12_spp[1,] == 2)

                        if f2_ == 0:
                            f2_ = 1
                        if f_2 == 0:
                            f_2 = 1

                        U1 = np.sum(D12_spp[0,]/n)
                        V1 = np.sum(D12_spp[1,]/m)
                        U3 = f_1/(2 * f_2)
                        V3 = f1_ / (2 * f2_)
                        U4 = np.sum( (D12_spp[0,]/n) * (D12_spp[1,] == 1 )  )
                        V4 = np.sum( (D12_spp[1,]/m) * (D12_spp[0,] == 1 )  )

                        if freq == "None":
                            U2 = (m-1)/m
                            V2 = (n-1)/n

                        if freq != "None":
                            if not isinstance(freq, np.ndarray):
                                print("Argument 'freq' must be a numpy ndarry")
                            if len(freq) != dat.shape[0]:
                                print("The length of freq must be same as sample #")

                            w = freq[i]
                            z = freq[j]
                            U2 = (z-1)/z
                            V2 = (w-1)/w

                        U = U1 + U2 * U3 * U4
                        V = V1 + V2 * V3 * V4
                        if U > 1:
                            U = 1.0
                        if V > 1:
                            V = 1.0

                    if method == "chao_jaccard":
                        sim0[j,i] = (U * V)/(U + V - U*V)
                    if method == "chao_sorensen":
                        sim0[j,i] = (2 * U * V)/(U+V)
                        print(sim0)

        sim = sim0[1,0]
        if sim == 1:
            sim = 0.99
    else:
        print("'{0}' is not valid similarity method, please choose one from {1}".format(method, avail_method))

    return(sim)


#
# a = [1,0,4,0,0,0,0,1]
# b = [2,1,3,0,0,1,0,6]
# make_matrix(a,b)
# similarity_trait(trait1=a, trait2 = b,method = "chao_sorensen")