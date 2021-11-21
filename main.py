# Test example for beesearch & BeeVolve

from BeeSearch import GoogleSearch
from BeeVolve import *

# Strand instanciation test

def usage():
    print("""
    Usage : 

    > python {} [-f|--file <file.tsp>] [-n|--niter <number>] [-p|--pmod <size>] [-r|--rate <number>] 
                {{-e|--execs <number>}} {{-i|--inittype [RGT|NNI]}} {{-s|--selection [RDM|BIN]]}} 
                {{-c|--crossover [UNI|ORD|USR|]}} {{-m|--mutation [SCR|INV]}} 
                {{-a|--alterpct [0-100]}} {{-x|--xovrpct [0-100]}} 
                {{-v|--verbose}} {{-w|--writestats}} {{-b|--bestonly}}

    Arguments 

    -f|--file               <file.tsp> : TSP instance file to be searched 

    -n|--niter              <integer>  : number of iterations before exiting 

    -p|--pmod               <integer>  : population size modulator

    -r|--rate               <number>   : probability for a mutation to occur


    Optional :

    -i|--inittype           [RGS|NNI] : Initialization mode, either RGS (Randomly Generated Strands) or 
                                        NNI (Nearest neighbor insertion) - default : RGT

    -s|--selection          [RDM|BIN|USR|FRD] : Type or selection performed from the matting pool, RDM (Random),
                                        BIN (BinaryTournament), USR (User-based), FRD (relation-based) - default : RDM

    -c|--crossover          [UNI|ORD|STG] : Crossover operator choice, either UNI (Uniform), ORD (Order-1) 
                                        of STG (Stigmer) crossover - default : UNI 
                                        UNI and ORD are automated test modes (fitness compared to original indexing)
                                        STG is production mode (fitness evaluated based on user decisions)
                                        Other mixed modes might be possible in the future.

                                        Note : ORD, STG not implemented yet ! 

    -m|--mutation           [SCR|INV] : Mutation operator choice, either SCR (Scrambling) or INV (Inversion)
                                        - default : SCR 
                                        
                                        Note : INV not implemented yet !

    -a|--alterpct           <integer> : determines in what proportion (percentage of genes)
                                        an individual will be mutated whenever a mutation occurs

    -x|--xovrpct            <integer> : determines the percentage of genes impacted by a crossover operation
                                        This will be used to compute the crossover section size 
                                        (default : 100% for Uniform)

    -v|--verbose            Display additional informations about the performance of the run and graph the fitness

    -w|--writestats         Write additional informations about the performance to a log and save the associated graph

    -b|--bestonly           Update an individual in the mating pool only if its replacement has better fitness

    """.format(sys.argv[0]))


def main(argv):
    if not os.path.exists("output"):
        os.mkdir("output")

    # Parameters with default values
    strandtype = "RGS"
    xovrtype = "SIM"
    slcttype = "BIN"
    mutntype = "SCR"
    bestonly = False
    strandsize = 21

    # Output parameters
    verbose = False
    outfile = None
    checkpointfile = None  # Unused yet, will save intermediate population state

    # parameters with unset checking
    infile = None
    psize = None
    niter = None
    rate = None
    apct = None
    xpct = None

    try:
        if len(argv) < 2:
            usage()
            sys.exit(2)
        opts, args = getopt.getopt(argv, "f:i:s:c:m:r:p:n:a:x:o:z:vb",
                                   ["file=", "inittype=", "selection=", "crossover=", "mutation=",
                                    "rate=", "psize=", "niter=", "alterpct=", "xovrpct=", "outfile=",
                                    "strandsize=", "verbose", "bestonly"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:

        try:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()
            elif opt in ("-r", "--rate"):
                rate = float(arg)
            elif opt in ("-p", "--psize"):
                psize = int(arg)
            elif opt in ("-n", "--niter"):
                niter = int(arg)
            elif opt in ("-z", "--strandsize"):
                strandsize = int(arg)
            elif opt in ("-a", "--alterpct"):
                if int(arg) in range(0, 100):
                    apct = int(arg)
                else:
                    print("[!] Invalid alteration percentage, will be set to None")
            elif opt in ("-x", "--xovrpct"):
                if int(arg) in range(0, 100):
                    xpct = int(arg)
                else:
                    print("[!] Invalid alteration percentage, will be set to None")
            elif opt in ("-v", "--verbose"):
                verbose = True

            elif opt in ("-w", "--outfile"):
                outfile = arg

            elif opt in ("-b", "--bestonly"):
                bestonly = arg

            elif opt in ("-f", "--file"):
                infile = arg

            elif opt in ("-i", "--inittype"):
                if arg in ("NNI", "RGS"):
                    strandtype = arg
                else:
                    print("[!] Unrecognized initialization type, defaulting to RGS - Randomly Generated Strands")
            elif opt in ("-s", "--selection"):
                if arg in ("RDM", "BIN"):  # TODO : add "STG", "REL"
                    slcttype = arg
                else:
                    print("[!] Unrecognized selection operator, defaulting to BIN - Binary Tournament")
            elif opt in ("-c", "--crossover"):
                if arg in ("SIM"): # TODO : add "STG"
                    xovrtype = arg
                else:
                    print("[!] Unrecognized crossover operator, defaulting to SIM - Similarity-Based Crossover")
            elif opt in ("-m", "--mutation"):
                if arg in ("SCR"): # TODO : add "INV"
                    mutntype = arg
                else:
                    print("[!] Unrecognized mutation operator, defaulting to SCR - Scrambling")

        except:
            print("[E] Some argument is invalid")

    if infile is None:
        print("[E] No TSP instance file provided. missing -f <file.tsp>")
        usage()
        sys.exit()
    if psize is None:
        print("[E] No population size provided. missing -p <size>")
        usage()
        sys.exit()
    if rate is None:
        print("[E] No mutation rate provided. missing -r <rate 0 to 1>")
        usage()
        sys.exit()
    if niter is None:
        print("[E] No number of iterations provided. missing -n <number>")
        usage()
        sys.exit()

    run_stat = []
    beeVolve = BeeVolve(_inputFile=infile,
                        _psize=psize,
                        _mutRate=rate,
                        _maxIter=niter,
                        _initType=strandtype,
                        _xovrType=xovrtype,
                        _slctType=slcttype,
                        _mutnType=mutntype,
                        _strandSize=strandsize,
                        _apct=apct,
                        _xpct=xpct,
                        _verbose=verbose,
                        _bestOnly=bestonly,
                        _outfile=outfile, _chkptFile=None)

    beeVolve.run()

    t_slct = 0
    t_xovr = 0
    t_mutn = 0
    for i in range(len(beeVolve.stat_board)):
        t_slct += beeVolve.stat_board[i][0]
        t_xovr += beeVolve.stat_board[i][1]
        t_mutn += beeVolve.stat_board[i][2]

    best_fitness = beeVolve.best.getStrandFitness()
    # TODO : Include implementation for user based fitness measurement
    run_stat.append([beeVolve.stat_inittime, t_slct, t_xovr, t_mutn, best_fitness])


if __name__ == "__main__":
    main(sys.argv[1:])


# Run example :
# python main.py -f test_urls.txt -n 200 -p 200 -r 0.01 -a 50 -z 21 -v


# Search instanciation test

# b = GoogleSearch("Economie bleue")
# results = b.getResults()
# for result in results:
#     print(result)

"""
(0, 'https://ec.europa.eu/commission/presscorner/detail/fr/ip_21_2341')
(1, 'https://www.un.org/africarenewal/fr/magazine/d%C3%A9cembre-2018-mars-2019/economie-bleue-une-opportunit%C3%A9-pour-l%E2%80%99afrique')
(2, 'https://fr.wikipedia.org/wiki/%C3%89conomie_bleue')
(3, 'https://www.banquemondiale.org/fr/topic/oceans-fisheries-and-coastal-economies')
[...]
(122, 'https://www.challenges.fr/')
(123, 'https://particulier.edf.fr/fr/accueil/gestion-contrat/options/ejp.html')
(124, 'https://www.linkedin.com/company/expertise-france')
"""