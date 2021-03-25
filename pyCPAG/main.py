#! /usr/bin/env python


from __future__ import generators

__author__ = "Liuyang Wang (wallacewly@gmail.com)"
__copyright__ = "Copyright (C) 2020 Liuyang Wang"
__license__ = "GPL 3.0"

soft_info = '''
Name:           CPAG.py
Purpose:        Cross Phenotype Aanalysis of GWAS
Author:         Liuyang Wang
Created on:     Thur May 25 2020
Copyright:      (c) Liuyang Wang
Licence:        GPL v 3.0
'''

import datetime, fnmatch
import sys
from _utils import *
#from heatmap import *

def gwasdb_cmp(args):
    """
    :param args: a list of parameters, type '-h' for more options
    :return: result, a pandas DataFrame
    """
    if len(sys.argv) == 2:
        sys.exit(0)

    gwasdb = query_gwasumdb()
    gwasdb.conndb()
    #gwasdb.search_traits(itrait="Interstitial lung disease", pcut=1e-5)

    all_subtypes = gwasdb.get_allsource();
    alltraits = gwasdb.get_alltrait()
    print("iCPAGdb has a total of {} GWAS traits/disease".format(len(alltraits)))

    # clin_gwastype = ["IBDgwas", 'NHGRI', 'Neale-B_UKBB_EUR_2017']  ## Here I define those 3 as clinical GWAS, all other are molecular traits
    clin_gwastype = ["IBDgwas", 'NHGRI']
    mole_gwastype = list(set(all_subtypes) - set(clin_gwastype))
    gwastypedict = {'mol_gwas': mole_gwastype, 'clin_gwas': clin_gwastype}

    subtypes = []; outdf = []; oname = ""

    nowdate = datetime.datetime.now()
    out_postsuffix = str(nowdate.year) + str(nowdate.month) + str(nowdate.day)

    if not args.subtype:
        print("please choose a subtype with --subtype")

    if args.pcut:
        pcut = float(args.pcut)
    else:
        pcut = 1e-5

    if len(args.subtype) == 1:
        assert isinstance(args.subtype, list)
        assert len(args.subtype) == 1

        if args.subtype[0] == 'mol_gwas':
            subtypes.extend(gwastypedict['mol_gwas'])
            print("CPAG analyzes molecular GWAS of " + ", ".join(gwastypedict['mol_gwas']))
        elif args.subtype[0] == 'clin_gwas':
            subtypes.extend(gwastypedict['clin_gwas'])
            print("CPAG analyzes clinical GWAS of " + ", ".join(gwastypedict['clin_gwas']))
        elif "," in args.subtype[0]:
            subtype_split = args.subtype[0].replace(" ", "").split(",")
            for itype in subtype_split:
                if itype not in all_subtypes:
                    sys.exit(str(itype) + " is not a valid GWAS subtype")
                subtypes.append(itype)
            print("CPAG analyzes the GWAS traits of {} ".format(str(args.subtype[0])))
        else:
            if args.subtype[0] not in all_subtypes:
                sys.exit(str(args.subtype[0]) + " is not a valid GWAS subtype")

            print("CPAG analyzes GWAS subtype: {} ".format(args.subtype[0]))
            subtypes.append(args.subtype[0])

        subtypeall = set(subtypes)
        subtypelist = []

        for itrait in subtypeall:
            if itrait == "H2P2":
                tmpdat = gwasdb.search_source(source=itrait, pcut=args.h2p2_pcut)
                # h2p2_prunedat = H2P2(pvcut = args.h2p2_pcut)
                # tmpdat = h2p2_prunedat.snpdict()
            elif itrait == "NHGRI":
                tmpdat = gwasdb.search_source(source=itrait, pcut=float(args.nhgri_pcut))
            else:
                tmpdat = gwasdb.search_source(source=itrait, pcut=pcut)

            subtypelist.append(tmpdat)

        snpdat = dict(ChainMap(*subtypelist))
        snpdat = {key: value for key, value in snpdat.items() if len(value) > 0}

        outdf = intra_trait(snpdat, ncpus=args.cpu,ldpop=args.ldpop, r2cut=args.ldr2,  cross_traits=False)
        oname = "-".join([x for x in subtypeall]).replace(" ", "").replace(",", "")
    else:
        if len(args.subtype) > 2:
            sys.exit("Please limit to compare two subtype. For more than 2 subtypes, use [--subtype A,B --subtype C] ")

        print("CPAG analyze GWAS traits between subtypes '{0}' and '{1}' ".format(args.subtype[0], args.subtype[1]))

        for i_subtype in args.subtype:
            if "," in i_subtype:
                subtype_split = i_subtype.replace(" ", "").split(",")
                for itype in subtype_split:
                    if itype == "mol_gwas" or itype == "clin_gwas":
                        sys.exit('Please run mol_gwas/clin_gwas separately, not use "," ')
                    if itype not in all_subtypes:
                        print(" {} is not a valid GWAS subtype".format(itype))
                        msg = 'select from {}, \nalternatively, select ["mol_gwas", "clin_gwas"] to choose all molecular GWAS traits or clinical GWAS traits '.format(
                            ", ".join(all_subtypes))
                        print(msg)
                        sys.exit(1)
                subtypes.append(subtype_split)
            elif i_subtype not in all_subtypes + ["mol_gwas", "clin_gwas"]:
                print(" {} is not a valid GWAS subtype".format(i_subtype))
                msg = 'Please select from {}, \nalternatively, select ["mol_gwas", "clin_gwas"] to choose all molecular GWAS traits or clinical GWAS traits '.format(
                    ", ".join(all_subtypes))
                print(msg)
                sys.exit(1)
            else:
                if i_subtype == "mol_gwas":
                    subtypes.append(gwastypedict['mol_gwas'])
                elif i_subtype == "clin_gwas":
                    subtypes.append(gwastypedict['clin_gwas'])
                else:
                    subtypes.append(i_subtype)

        subtypeall = subtypes
        subtypelist = []

        for itrait in subtypeall:
            if isinstance(itrait, list):
                tmp_sublist = []

                for sub_itrait in itrait:
                    if sub_itrait == "H2P2":
                        tmpdat = gwasdb.search_source(source=sub_itrait, pcut=args.h2p2_pcut)
                        # h2p2_prunedat = H2P2(pvcut = args.h2p2_pcut)
                        # tmpdat = h2p2_prunedat.snpdict()
                    elif sub_itrait == "NHGRI":
                        tmpdat = gwasdb.search_source(source=sub_itrait, pcut= float(args.nhgri_pcut))
                    else:
                        tmpdat = gwasdb.search_source(source=sub_itrait, pcut=pcut)
                    tmp_sublist.append(tmpdat)
                subtypelist.append(dict(ChainMap(*tmp_sublist)))
            else:
                if itrait == "H2P2":
                    tmpdat = gwasdb.search_source(source=itrait, pcut=args.h2p2_pcut)
                    # h2p2_prunedat = H2P2(pvcut = args.h2p2_pcut)
                    # tmpdat = h2p2_prunedat.snpdict()
                elif itrait == "NHGRI":
                    tmpdat = gwasdb.search_source(source=itrait, pcut=float(args.nhgri_pcut))
                else:
                    tmpdat = gwasdb.search_source(source=itrait, pcut=pcut)

                subtypelist.append(tmpdat)

        snpdat = subtypelist
        outdf = intra_trait(snpdat, ncpus=args.cpu,ldpop=args.ldpop, r2cut=args.r2cut,  cross_traits=True)

        oname_merge = []
        for x in subtypeall:
            if isinstance(x, list):
                oname_merge.append("-".join(x))
            else:
                oname_merge.append(x)

        oname = "_vs_".join([x for x in oname_merge]).replace(" ", "").replace(",", "")

    if "H2P2" in subtypes:
        if len(subtypes) == 1:
            out_pfx = "./output/" + "cpag_output_" + "H2P2pcut" + str(args.h2p2_pcut) + "_" + out_postsuffix + ".out.csv"
        else:
            out_pfx = "./output/" + "cpag_output_" + "H2P2pcut" + str(args.h2p2_pcut) + "-" + oname + "pcut" + str(
                pcut) + "_" + out_postsuffix + ".out.csv"
    else:
        out_pfx = "./output/" + "cpag_output_" + oname + "_pcut" + str(pcut) + "_" + out_postsuffix + ".out.csv"

    if outdf.shape[0] < 1:
        sys.exit("There is no overlapping SNPs! Finished!")

    if args.outfile:
        outdf.to_csv(args.outfile, index=False)
        print(f"Result is saved into {args.outfile}.")
    else:
        outdf.to_csv(out_pfx, index=False)
        print(f"Result is saved into {out_pfx} .")


def load_dat(files, args):
    """
    :param indir: str. infile dir.
    :param files: str. infiles
    :param delimiter: str. '\t' for tab; ',' for csv; ' ' or '\s' for white space
    :param SNPcol: str. Column name for SNP rsID
    :param Pcol: str. Column name for P value
    :param pcut: float. Exclude SNPs with p value > pcut.
    :return: dict. Dict keys = file name, and SNP list are values
    """

    dat_dict = {};
    ldproxy_snps = {}

    if len(files) > 1:
        for i_file in files:
            if i_file.endswith("gz"):
                idat_chunk = pd.read_csv(i_file, sep=args.delimitor,  compression='gzip', iterator=True,
                                         chunksize=100000)
            else:
                idat_chunk = pd.read_csv(i_file, sep=args.delimitor, iterator=True, chunksize=100000)
            idat1 = pd.concat([chunk[chunk[args.pcol] < float(args.usrpcut)] for chunk in idat_chunk])
            idat1.fillna("NA", inplace=True)

            if args.ldclump:
                idat, ldproxy_snps = ld_clump(idat1, args)
                dat_dict[i_file] = idat[args.snpcol].values.tolist()
            else:
                dat_dict[i_file] = idat1[args.snpcol].values.tolist()
    else:
        if files[0].endswith("gz"):
            idat_chunk = pd.read_csv(files[0], sep=args.delimitor, compression='gzip',  iterator=True,
                                     chunksize=100000)
        else:
            idat_chunk = pd.read_csv(files[0], sep=args.delimitor,iterator=True, chunksize=100000)

        usr_dat1 = pd.concat([chunk[chunk[args.pcol] <= float(args.usrpcut)] for chunk in idat_chunk])
        # usr_dat1 = pd.read_csv(files[0], sep=args.delimitor)
        usr_dat1.fillna("NA", inplace=True)

        if args.ldclump:
            usr_dat, ldproxy_snps = ld_clump(usr_dat1, args)
            dat_dict[files[0]] = usr_dat[args.snpcol].values.tolist()
        else:
            dat_dict[files[0]] = usr_dat1[args.snpcol].values.tolist()

    if len(dat_dict) < 1:
        sys.exit("No SNPs were kept after filtering P values of " + str(args.usrpcut))
    # else:
    #     snplens = [len(v) for k, v in dat_dict.items()]
    #     # snplen = sum(snplens)
    #     # print(str(snplen) + " SNPs are left for cross-phenotype analysis ")

    return ([dat_dict, ldproxy_snps])


def user_gwas_cmp(args):
    """
    allow user to upload own GWAS dataset, and compute the similarity against iCPAGdb
    :param args:
    :return:
    """

    if len(sys.argv) == 2:
        sys.exit(0)

    gwasdb = query_gwasumdb()
    gwasdb.conndb()

    all_subtypes = gwasdb.get_allsource()
    alltraits = gwasdb.get_alltrait()
    print("iCPAGdb has a total of {} GWAS traits/disease".format(len(alltraits)))
    print(f"Current GWAS data source in CPAG database are: {all_subtypes} ")

    clin_gwastype = ["IBDgwas", 'NHGRI', 'Neale-B_UKBB_EUR_2017']  ## Here I define those 3 as clinical GWAS, all other are molecular traits
    mole_gwastype = list(set(all_subtypes) - set(clin_gwastype))
    gwastypedict = {'mol_gwas': mole_gwastype, 'clin_gwas': clin_gwastype}

    file_dir = os.getcwd() if args.indir == None else args.indir

    nowdate = datetime.datetime.now()
    out_postsuffix = str(nowdate.year) + str(nowdate.month) + str(nowdate.day)

    if not os.path.exists(file_dir):
        sys.exit("Folder {} does not exist" % file_dir)

    all_infiles = []

    if not args.infile:
        sys.exit("Please choose a input file with --infile")

    if "," in args.infile:
        my_infiles = args.infile.replace(" ", "").split(",")
        for ifile in my_infiles:
            fileexist = check_file(file_dir, ifile)
            if fileexist:
                all_infiles.append(os.path.join(file_dir, ifile))

    elif "*" in args.infile:
        for file in os.listdir(args.indir):
            if fnmatch.fnmatch(file, args.infile):
                fileexist = check_file(file_dir, file)
                if fileexist:
                    all_infiles.append(os.path.join(file_dir, file))
    else:
        filexist = check_file(file_dir, args.infile)
        if not filexist:
            sys.exit(" {} file does not exist.".format(str(args.infile)))
            # all_infiles.append(os.path.join(file_dir,args.infile))
        all_infiles.append(args.infile)

    user_dat, usr_proxysnpdict = load_dat(all_infiles, args)

    pcut = args.pcut
    subtypeall = []
    if args.subtype:
        if args.subtype not in all_subtypes + ['mol_gwas','clin_gwas','all','All','ALL']:
            sys.exit(str(args.subtype) + " is not a valid GWAS subtype")
        else:
            print("Loading {} GWAS traits in CPAG db ...".format(args.subtype))
            if args.subtype in ["all","All", "ALL"]:
                subtypeall = all_subtypes
            elif args.subtype == "mol_gwas":
                subtypeall.extend(gwastypedict['mol_gwas'])
            elif args.subtype == "clin_gwas":
                subtypeall.extend(gwastypedict['clin_gwas'])
            else:
                subtypeall.extend([args.subtype])

    subtypelist = []
    for itrait in subtypeall:
        if itrait == "H2P2":
            tmpdat = gwasdb.search_source(source=itrait, pcut = args.h2p2_pcut)
            # h2p2_prunedat = H2P2(pvcut=args.h2p2_pcut)
            # tmpdat = h2p2_prunedat.snpdict()
        elif itrait == "NHGRI":
            tmpdat = gwasdb.search_source(source=itrait, pcut=float(args.nhgri_pcut))
        else:
            tmpdat = gwasdb.search_source(source=itrait, pcut = pcut)

        subtypelist.append(tmpdat)

    snpdat_db = dict(ChainMap(*subtypelist))
    snpdat_db = {k: v for k, v in snpdat_db.items() if len(v) > 0}
    print("A total of {} traits/disease were load from CPAG db v2.0".format(len(snpdat_db)))

    ## No need for pre-filtering as the speed is good
    # if usr_proxysnpdict:
    #     # pre-filter GWAS database
    #     ## filter built-in GWAS as if there is no SNPs in a GWAS overlapping with usr's proxy + leading SNPs, this traits can be removed
    #     usr_proxysnpset = set(usr_proxysnpdict["leadingSNPs"]).union(set(usr_proxysnpdict["proxySNPs"]))
    #     del_trait = []
    #     for k, v in snpdat_db.items():
    #         xhit = set(v).intersection(usr_proxysnpset)
    #         if len(xhit) < 1:
    #             del_trait.append(k)
    #     if len(del_trait) > 0:
    #         for i_deltrait in del_trait:
    #             # del snpdat_db[i_deltrait]
    #             continue
    # print("after ld clumpping", len(snpdat_db))

    snpdat = [user_dat,snpdat_db]
    outdf = intra_trait(snpdat, ncpus=args.cpu,ldpop=args.ldpop, r2cut=args.ldr2, cross_traits=True)

    outdf['Trait1'] = str(args.usr_phename)

    if not outdf.empty:
        out_pfx = "./output/" + "cpag_output_usrGWASpcut" + str(args.usrpcut) + "_cpagDBpcut" + str(pcut) + "_" + out_postsuffix + ".out.csv"
        if args.outfile:
            outdf.to_csv(args.outfile, index=False)
            print("Result is saved into {}.!".format(args.outfile))
        else:
            outdf.to_csv(out_pfx, index=False)
            print("Result is saved into {}.!".format(out_pfx))

def post_analysis(args):

    if args.ontology:
        if not args.infile:
            print("Please choose the input file with '--infile' ")
            sys.exit(0)

        import anno_parent_efo

        annocol = args.annocols
        print(f"\nAdding EBI EFO parent ontology term for {annocol} ")

        dat = pd.read_csv(args.infile, na_filter="NA")
        if not args.outfile:
            outfile = args.infile.replace(".csv","_addOntology.csv")
        else:
            outfile =args.outfile

        annocols = annocol.split(",")

        datfin = pd.DataFrame()

        if len(annocols) > 1:
            dtmp = [dat]
            for i in range(len(annocols)):
                if annocols[i] not in dat.columns:
                    print(f"Column {annocols[i]} is not a valid column in {args.infile}")
                    sys.exit(0)

                rout = anno_parent_efo.add_oto(idat=dat, annocol=annocols[i])
                r_df = pd.DataFrame(rout, columns =[annocols[i] + "_EFO",annocols[i] + "_ParentTerm",annocols[i] + "_ParentEFO"])
                dtmp.append(r_df)
            datfin = pd.concat(dtmp, axis = 1) # ignore_index=True,
        else:
            if annocols[0] not in dat.columns:
                print(f"Column {annocols[0]} is not a valid column for {args.infile}")
                sys.exit(0)
            rout = anno_parent_efo.add_oto(idat=dat, annocol=annocols[0])
            r_df = pd.DataFrame(rout, columns =[annocols[0] + "_EFO",annocols[0] + "_ParentTerm",annocols[0] + "_ParentEFO"])
            datfin = pd.concat([dat, r_df], axis = 1)

        datfin.to_csv(outfile, index=False)

    if args.heatmap:
        print("Not implement this function yet")
        #heatmapplot(args.infile)

def demo(args):
    """
    demo of running example
    """
    if args.rungwas:
        print("Here are available options:")
        print("--> running within a subtype")
        run_cmd = "python cpag.py cpagdb --subtype 'H2P2' \n "
        print(run_cmd)

        print("--> running within two subtype")
        run_cmd = "python cpag.py cpagdb --subtype 'H2P2, IBDgwas' \n  "
        print(run_cmd)

        print("--> running within multiple subtype")
        run_cmd = "python cpag.py cpagdb --subtype 'H2P2, IBDgwas,BloodCytokine' \n  "
        print(run_cmd)

        print("--> running within all molecular GWAS traits")
        run_cmd = "python cpag.py cpagdb --subtype 'mol_gwas' \n "
        print(run_cmd)

        print("--> running within all clinical GWAS traits")
        run_cmd = "python cpag.py cpagdb --subtype 'clin_gwas' \n  "
        print(run_cmd)

        print("--> running between two subtypes")
        run_cmd = "python cpag.py cpagdb --subtype 'H2P2' --subtype 'IBDgwas'  \n "
        print(run_cmd)

    if args.usrgwas:
        print("--> running usr's own GWAS dataset against H2P2")
        run_cmd = "python cpag.py usr_gwas --infile myown_gwas.txt --query 'H2P2' \n "
        print(run_cmd)

        print("--> running usr's own GWAS dataset against all GWAS traits in iCPAGdb")
        run_cmd = "python cpag.py usr_gwas --infile myown_gwas.txt --Pcol 'Pvalue' --SNPcol 'SNPs' --query 'all' \n "
        print(run_cmd)

        print("--> running usr's own GWAS dataset against molecular GWAS traits")
        run_cmd = "python cpag.py usr_gwas --infile myown_gwas.txt --Pcol 'Pvalue' --SNPcol 'SNPs' --query 'mol_gwas' \n "
        print(run_cmd)

        print("--> running usr's own GWAS dataset against clinical GWAS traits")
        run_cmd = "python cpag.py usr_gwas --infile myown_gwas.txt --Pcol 'Pvalue' --SNPcol 'SNPs' --query 'mol_gwas' \n "
        print(run_cmd)


def main(prog=None):
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Cross-Phenotype Analysis of GWAS", add_help=True,
                                     prog=prog if prog else sys.argv[0])
    root_parser = parser.add_subparsers()

    inter_gwas = root_parser.add_parser('cpagdb', help='Cross-phenotype analysis of GWAS traits in CPAG database')
    inter_gwas.add_argument('-T', '--subtype', dest="subtype", action='append',
                            # choices=["H2P2","BloodCytokine","BloodMetabolites", "BloodXenobiotic","IBDgwas","UrineMetabolites", "mol_gwas", "clin_gwas"],
                            help='use [--subtype A] to analyze traits within subtype A;  or '
                                 '[--subtype "A,B"] for analyzing GWAS traits within both subtypes [A + B]; or '
                                 '[--subtype "A" --subtype "B"] to analyze GWAS traits in A against B')
    inter_gwas.add_argument('--Pcut', type=float, dest="pcut", action="store", default=5e-8, help="filtering out SNPs")
    inter_gwas.add_argument('--NHGRI-Pcut', type=float, dest="nhgri_pcut", action="store", default=5e-8, help="filtering out SNPs")
    inter_gwas.add_argument('--H2P2-Pcut', type=float, dest="h2p2_pcut", action="store", default=1e-5,
                            help="filtering out SNPs for H2P2")
    inter_gwas.add_argument('--lddb-pop', type=str, dest="ldpop", default='EUR',
                            help="LD reference population, default: 'EUR', available for ['EUR', 'AFR', 'EAS'] ")
    inter_gwas.add_argument('--lddb-r2', type=float, dest="ldr2", default=0.4,
                           help="LD database, default: '0.4', available for [0.2, 0.4, 0.8] ")
    inter_gwas.add_argument('--threads', type=int, dest="cpu", action="store", default=1,
                            help="Set multiple CPUs/threads. Default: 1")
    inter_gwas.add_argument('-O', '--outfile', type=str, dest="outfile", action="store", default=None,
                            help="output csv file name")
    inter_gwas.set_defaults(func=gwasdb_cmp)

    user_gwas = root_parser.add_parser('usr-gwas',
                                       help='Compare user GWAS against CPAGdb traits. Type -h for more option')
    user_gwas.add_argument('--indir', type=str, dest="indir", default=None,
                           help="input files directory. Default: current folder")
    user_gwas.add_argument('--infile', type=str, dest="infile",
                           help="input files, use '--infile A_file.txt, B_file.txt' or "
                                "'--infile *_file.txt' for multiple files")
    user_gwas.add_argument('--delimitor', dest="delimitor", default='\t',
                           help="Delimiter to use. e.g. '\t' for tab, ',' for csv, ' ' for white space ")
    user_gwas.add_argument('--SNPcol', type=str, dest="snpcol", default='SNP', help="SNP rsID column name in the file")
    user_gwas.add_argument('--Pcol', type=str, dest="pcol", default='P',
                           help="Column name for p value")  # action='store_true')
    user_gwas.add_argument('--usr-pcut', type=float, dest="usrpcut", default=1e-5,
                           help="p cutoff for user's GWAS. Default: 1e-5")
    user_gwas.add_argument('--NHGRI-Pcut', type=float, dest="nhgri_pcut", action="store", default=5e-8,
                           help="filtering out SNPs")
    user_gwas.add_argument('--usr-pheno-name', type=str, dest="usr_phename", default="User_trait",
                           help="User phenotype name for Trait1 in outfile. Default: User_trait")
    user_gwas.add_argument('--H2P2-Pcut', type=float, dest="h2p2_pcut", action="store", default=1e-5,
                           help="filtering out SNPs for H2P2. Default: 1e-5")
    user_gwas.add_argument('--cpagdb-pcut', type=float, dest="pcut", default=5e-8,
                           help="p cutoff for CPAG db GWAS. Default: 5e-8")
    user_gwas.add_argument('--ld-clump', type=float, dest="ldclump", default=1,
                           help="perform LD clumping, choose [1/0]. Default: 1")
    user_gwas.add_argument('--lddb-pop', type=str, dest="ldpop", default="EUR",
                           help="Population used for LD clumping, default: 'EUR', available for ['EUR', 'AFR', 'EAS'] ")
    user_gwas.add_argument('--lddb-r2', type=float, dest="ldr2", default=0.4,
                           help="LD database, default: '0.4', available for [0.2, 0.4, 0.8] ")
    user_gwas.add_argument('--ld-clump-p1', type=float, dest="ldclump_p1", default=1e-5,
                           help="parameters for LD clumping p1, default: 1e-5")
    user_gwas.add_argument('--ld-clump-p2', type=float, dest="ldclump_p2", default=1,
                           help="parameters for LD clumping p2, default: 1")
    user_gwas.add_argument('--ld-clump-r2', type=float, dest="ldclump_r2", default=0.4,
                           help="LD clump r-squared cutoff, default: 0.4")
    user_gwas.add_argument('--ld-clump-kb', type=float, dest="ldclump_kb", default=5000,
                           help="LD clump window size, default: 5000")
    user_gwas.add_argument('--threads', type=int, dest="cpu", action="store", default=1,
                           help="Set multiple CPUs/threads. Default: 1")
    user_gwas.add_argument('--querydb', dest="subtype", default="all",
                           # choices=["H2P2","BloodCytokine","BloodMetabolites", "BloodXenobiotic","IBDgwas","UrineMetabolites", "mol_gwas", "clin_gwas"],
                           help='query single GWAS source, or all traits. Default: all')
    user_gwas.add_argument('-O', '--outfile', type=str, dest="outfile", action="store", default=None,
                           help="output csv file name")
    user_gwas.set_defaults(func=user_gwas_cmp)

    post_parser = root_parser.add_parser('post_analysis', help='post analysis for heatmap/ontology. Not full developed yet')
    post_parser.add_argument("--heatmap", dest='heatmap', action='store_true', help = "plot heatmap")
    post_parser.add_argument("--oneline-heatmap", dest='oneheatmap', action='store_true', help = "plot oneline heatmap")
    post_parser.add_argument("--network", dest='network', action='store_true', help = "plot network")
    post_parser.add_argument("--anno-ontology", dest='ontology', action='store_true', help="add EBI ontology EFO to disease/trait")
    post_parser.add_argument("--anno-cols", dest='annocols', default="Trait2", help="columns to annotate, default 'Trait2', or with multiple columns 'Trait1,Trait2' ")
    post_parser.add_argument("--infile", dest='infile', type=str, default=None,
                           help="A csv file. Either can be the icpagdb outfile, or any csv file has 'Trait1' and 'Trait2' column")
    post_parser.add_argument('-O','--outfile',  type=str, dest='outfile', default="icpagdb_post",
                             help="output filename")
    post_parser.set_defaults(func=post_analysis)

    demo_parser = root_parser.add_parser('demo', help='Demo/example')
    demo_parser.add_argument("--run-cpagdb", dest='rungwas', action='store_true')
    demo_parser.add_argument("--run-usrgwas", dest="usrgwas", action='store_true')
    demo_parser.set_defaults(func=demo)

    parser.add_argument('--version', action='version', version='%(prog)s 2.0')

    if len(sys.argv) == 1:
        print(parser.format_help())
        sys.exit(0)

    if len(sys.argv) == 2:
        """ help function for sub-function """
        subparsers_mod = [module for module in parser._actions if isinstance(module, argparse._SubParsersAction)]
        for i_subparsers in subparsers_mod:
            for choice, subparser in i_subparsers.choices.items():
                if sys.argv[1] == choice:
                    print(subparser.format_help())
                    sys.exit(0)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    start = datetime.datetime.now()
    main(prog="Cross Phenotype Analysis of GWAS")
    end = datetime.datetime.now()
    time_elapse = end - start
    print('Analysis is done! Total elapsed time (h:m:s):', time_elapse)
