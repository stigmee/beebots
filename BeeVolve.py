###############################################################################
## Stigmee: A 3D browser and decentralized social network.
## Copyright 2021 Duron Alain <duron.alain@gmail.com>
##
## This file is part of Stigmee.
##
## Project : Stigmee BeeBot
## Version : 0.0-1
## Date : 20-11-2021
## Author : Alain Duron
## File : BeeVolve.py
##
## Stigmee is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

# TODO:
#   Add the required structures to handle additional metadata within a Strand
#     (tags, pertivirality score, creator's hash adress...)
#   Add documentation to the header of each function
#   Add methods (or better, a separate module ?) to handle the database calls
#   Add separated execution modes : test (web scrapping) and run (using stigmergic functions)
#   Decouple the strand size (number of urls shown to user) and the total evaluation size (all urls)

import random
import sys
import os
import getopt
from BeeStrand import *
import timeit
import operator
import matplotlib.pyplot as mpl

from datetime import datetime

"""
At initialization, BeeVolve receives the result of a search, or a collection of urls from a database.
or static file. this collection is succeptible to contain a large number of url and therefore 
should be split into subsequent strands to be curated by users.

Splitting the data into several partitions containing different URLs is NOT a viable solution, as I experienced... 
So the evaluation could be done over the entire url set (for example 125 URls) even if the strand only shows 
20 of them to the user. 
"""

class BeeVolve:

    def __init__(self, _inputFile, _psize, _mutRate, _maxIter,
                 _initType, _xovrType, _slctType, _mutnType, _strandSize=21,
                 _apct=None, _xpct=None, _verbose=False, _bestOnly=False,
                 _outfile=None, _chkptFile=None):
        """
        :param _inputFile:
        :param _psize:
        :param _mutRate:
        :param _maxIter:
        :param _initType:
        :param _xovrType:
        :param _slctType:
        :param _mutnType:
        :param _strandSize:
        :param _apct:
        :param _xpct:
        :param _verbose:
        :param _bestOnly:
        :param _outfile:
        :param _chkptFile:
        """

        self.verbose = _verbose         # Wether messages should be displayed or not
        self.outfile = _outfile         # Wether a log should be kept
        self.inputFile = _inputFile     # If a checkpoint file is provided
        self.chkptFile = _chkptFile     # Checkpoint file used to save the actual data state (cache)
        idate = datetime.now()
        self.idate_str = idate.strftime("%Y%m%d_%H%M%S_")

        """
        Parameters and general variables
        """

        self.population = []  # The population of strand instances
        self.matingPool = []  # The pool used for strand pairing / evaluation
        self.best = None  # The fitest strand yet
        self.popSize = _psize  # The population size
        self.strandSize = _strandSize  # The number of URLs to load
        self.mutationRate = _mutRate  # Mutation rate (statistical)
        self.maxIter = _maxIter  # Maximum number of iteration (Testing)
        self.iteration = 0
        self.data = {}

        # Algorythm parameters
        self.initType = _initType   # RGS (Randomly generated Strand) or NNI (Nearest Neighbour)
        self.xovrType = _xovrType   # SIM (Similarity-based uniform) or ORD (order-1)
        self.slctType = _slctType   # RDM (Random) or BIN (Binary selection tournament)
        self.mutnType = _mutnType   # SCR (Scramble) or INV (Inversion
        self.apct = _apct           # The percentage of genes to consider when mutation occurs
        self.xpct = _xpct           # The section size to consider when crossover occurs
        # Allow replacement only when population individual is better than the one in the previous pool
        self.bestonly = _bestOnly

        # performance data
        self.stat_inittime = 0
        self.stat_board = []

        # Reads data file
        self.readUrlList()

        # output management
        self.iprint("[c] Num of iterations set to : {}".format(self.maxIter))
        self.iprint("[c] Population size set to : {}".format(self.popSize))
        self.iprint("[c] Mutation rate set to : {}".format(self.mutationRate))
        self.iprint("[c] Initialization type defined to : {}".format(self.initType))
        self.iprint("[c] Selection operator defined to : {}".format(self.slctType))
        self.iprint("[c] Crossover operator defined to : {}".format(self.xovrType))
        if self.xpct is not None:
            self.iprint("[c]     Crossover section size set to : {}".format(self.xpct))
        self.iprint("[c] Mutation operator defined to : {}".format(self.mutnType))
        self.iprint("[c]     Mutation Alter percentage set to : {}".format(self.apct))

        # Initialize and keep timed statistics
        beginInit = timeit.default_timer()
        self.initPopulation()
        endInit = timeit.default_timer()
        self.stat_inittime = endInit - beginInit
        self.iprint("[s] init time : {}".format(self.stat_inittime))

    def __del__(self):

        if self.outfile is not None: self.outfile.close()

    def iprint(self, message):
        """
        Print destination and display depending on output parameters
        """
        if self.outfile is not None:
            print(message, file=self.outfile)
        if self.verbose:
            print(message)

    def readUrlList(self):
        """
        Reading an list of URLs to be split into multiple strands for pupulation initialization
        """
        file = open(self.inputFile, 'r')
        with open(self.inputFile) as file:
            strandlines = [next(file) for l in range(self.strandSize)]
        file.close()
        self.data = {}
        for line in strandlines:
            (id, url) = line.split()
            self.data[int(id)] = url


    def initPopulation(self):
        """
        Initialize the population of strands for the algorithm
        """

        # 1. initialize the population up to popSize
        for j in range(0, self.popSize):

            strand = Strand(_strandSize=self.strandSize,
                            _origGenes=self.data,
                            _childGenes=None, _nni=False)
            strand.computeStrandFitness()
            self.population.append(strand)

        # Determine the initial population size
        self.popSize = len(self.population)
        self.iprint("Initial population size : {}".format(self.popSize))

        # Determine the best initial strand. In test mode, this uses index distance
        # and inversion counts only as no user is involved in the process.
        self.best = self.population[0].copy()
        for s in self.population:
            if self.best.getStrandFitness() > s.getStrandFitness():
                self.best = s.copy()

        self.iprint("Best initial sorting: {}".format(self.best.getStrandFitness()))

    def updateBest(self, candidate):
        if self.best is None or candidate.getStrandFitness() < self.best.getStrandFitness():
            self.best = candidate.copy()
            self.iprint("iteration: {} - best: {}".format(self.iteration,
                                                          self.best.getStrandFitness()))

    def stigmerSelection(self, sA, sB):
        # TODO:
        #   In the selection mode, the user picks up a predefined set of n strands
        #   to be mixed up and evaluated together
        #   evaluating mode than 2 strands together might reduce the curating process time
        pass

    def RelationalSelection(self, sA, sB):
        # TODO:
        #   For this mode, a stigmer get strands assigned based on relational information
        #   This can only be a valuable mode when the stigmer network have accumulated enough
        #   relational information
        pass

    def randomSelection(self):
        """
        Random (uniform) selection of two strands
        """
        sA = self.matingPool[random.randint(0, self.popSize - 1)]
        sB = self.matingPool[random.randint(0, self.popSize - 1)]
        return [sA, sB]

    def binaryTournamentSelection(self):
        """
        This method randomly selects 2 strands twice and match them,
        returning the fitest of each pair
        """
        winingStrands = []
        for _ in range(2):
            sA = self.matingPool[random.randint(0, self.popSize - 1)]
            sB = self.matingPool[random.randint(0, self.popSize - 1)]
            # Note : for strands, test evaluation is a minimization problem
            if sA.getStrandFitness() < sB.getStrandFitness():
                winingStrands.append(sA)
            else:
                winingStrands.append(sB)
        return winingStrands

    def stigmerCrossover(self, sA, sB):
        # TODO:
        #   For this mode, we don't care about actual sorting performance over the original index
        #   We only care about which one is rejected by the user based on the proposed url sets.
        #   Urls should be mixed up and proposed all at once, so that user is not influenced when selecting
        #   Useful work of curating the network can be rewarded
        #   It could also be interesting to check if the performance is close to the original
        #   sort (indicative of whether the data source is good at proposing good sets or urls)
        #   To be discussed with the dev because the way to implement this method impacts the interface
        pass

    def similarityBasedCrossover(self, sA, sB):
        """
        For this specific algorithm, strands can be initialized from different urls and therefore
        a uniform crossover is not possible over the whole strand. So instead we are going to
        implement a variant that i call similarity-based crossover

        Edit : for now i've restricted the number of Urls to the strand size

        1. I select only the similar entries accross both parent strandSize and create a sub index
        2. I perform a full uniform crossover on the sub index and create the childs accordingly

        When applied to a population where all individuals share the same genes, it is basically a uniform crossover.
        """

        commonset = set(sA.strandGenes) & set(sB.strandGenes)
        # common elements as they appear in strand A
        sAsub = sorted(commonset, key=lambda x: sA.strandGenes.index(x))
        # common elements as they appear in strand B
        sBsub = sorted(commonset, key=lambda x: sB.strandGenes.index(x))

        # If xpct is set, we define a section start and a section size of xpct
        # We set only 1s outsize of the section
        if self.xpct is not None:
            template = [1 for i in range(len(commonset))]
            sectionsize = int(len(commonset) * self.xpct / 100)
            sectionstart = random.randint(0, len(commonset) - sectionsize - 1)
            # set the section
            for i in range(sectionstart, sectionstart + sectionsize):
                template[i] = random.randint(0, 1)
        else:
            template = [random.randint(0, 1) for i in range(len(commonset))]

        remains = []
        subchildGenes = [0 for i in range(len(commonset))]

        # append from parent 1 when template[i] == 1
        for i in range(0, len(commonset)):
            if template[i] == 1:
                subchildGenes[i] = sAsub[i]
            else:
                remains.append(sAsub[i])

        # create sorted list of remaining items in sA ordered by their appearance in sB
        ordered_remains = [gA for gB in sBsub for gA in remains if gA == gB]

        pos = 0
        for i in range(0, len(commonset)):
            if template[i] == 0:
                subchildGenes[i] = ordered_remains[pos]
                pos += 1

        # Now that we got the subchild genes list, we reassemble them in either childA or childB
        childGenes1 = sA.strandGenes
        childGenes2 = sB.strandGenes
        for i in range(0, len(commonset)):
            childGenes1[sA.strandGenes.index(sAsub[i])] = subchildGenes[i]
            childGenes2[sB.strandGenes.index(sBsub[i])] = subchildGenes[i]

        # 2 possible child generation
        # Get the part of data associated with the child
        childdata1 = {}
        childdata2 = {}
        for i in range(0, len(childGenes1) ):
            childdata1[childGenes1[i]] = self.data[i]
        for i in range(0, len(childGenes2) ):
            childdata2[childGenes2[i]] = self.data[i]

        # Now we got a final child to propagate
        child1 = Strand(self.strandSize, childdata1, childGenes1)
        child2 = Strand(self.strandSize, childdata2, childGenes2)

        child1.computeStrandFitness()
        child2.computeStrandFitness()

        if child1.getStrandFitness() < child2.getStrandFitness():
            self.updateBest(child1)
            return child1
        else:
            self.updateBest(child2)
            return child2

    def scrambleMutation(self, strand):
        """
        This will select a subset of strandGenes and randomly shuffle them
        """
        if random.random() > self.mutationRate:
            return

        # Define a mutation section size
        if self.apct is not None:
            sectionsize = int(self.strandSize * self.apct / 100)
        else:
            sectionsize = random.randint(1, self.strandSize - 1)

        # Define a mutation section start
        sectionstart = random.randint(0, self.strandSize - sectionsize)

        # mutate by shuffling all strandGenes of the section
        subset = strand.strandGenes[sectionstart:sectionstart + sectionsize]
        random.shuffle(subset)
        strand.strandGenes[sectionstart:sectionstart + sectionsize] = subset

        strand.computeStrandFitness()
        self.updateBest(strand)

    def poolFitness(self):
        """
        get the best fitness of the mating pool.
        """
        bestStrandFitness = self.matingPool[0].getStrandFitness()
        for i in range(1, self.popSize):
            fit = self.matingPool[i].getStrandFitness()
            if fit < bestStrandFitness: bestStrandFitness = fit
        return bestStrandFitness

    def updateMatingPool(self):
        """
        Updating the mating pool before creating a new generation
        """
        # In bestOnly mode, an item only gets replaced in the matting pool if the replacement is fitter
        # This allows the pool to converge much faster
        if self.bestonly == False:
            self.matingPool = []
            for strand in self.population:
                self.matingPool.append(strand.copy())
        else:
            if len(self.matingPool) == 0: self.matingPool = []
            for i in range(0, self.popSize):  # self.population:
                if len(self.matingPool) == self.popSize:
                    if self.population[i].getStrandFitness() < self.matingPool[i].getStrandFitness():
                        self.matingPool[i] = self.population[i].copy()
                else:
                    self.matingPool.append(self.population[i].copy())

    def newGeneration(self):
        """
        Creating a new generation
        1. Selection
        2. Crossover
        3. Mutation
        """
        stat_gen = [0, 0, 0, 0, 0]

        for i in range(0, len(self.population)):

            # Depending on the slctType selected :
            s1 = timeit.default_timer()
            if self.slctType == "RDM":
                parent1, parent2 = self.randomSelection()
            else:  # "BIN":
                parent1, parent2 = self.binaryTournamentSelection()
            # TODO:
            #    add self.StigmerSelection() : User pick up of instances ("USR")
            #    add self.RelationalSelection() : Depending of friends preference on same topic ("FRD")

            # Depending on the xovrType selected :
            s2 = timeit.default_timer()
            child = self.similarityBasedCrossover(parent1, parent2)
            # TODO:
            #    add self.stigmerCrossover() crossover implementation (user-based)

            # Depending on the mutnType selected :
            s3 = timeit.default_timer()
            self.scrambleMutation(child)
            # TODO:
            #    add self.inversionMutation() type

            s4 = timeit.default_timer()

            # Updating the population with the child
            # We replace the strand that lost the round
            self.population[i] = child

            # filling in stats
            stat_i = [s2 - s1, s3 - s2, s4 - s3, 0]
            stat_gen = list(map(operator.add, stat_gen, stat_i))

        # append best fitness so far to the statistics
        stat_gen[3] = self.best.getStrandFitness()
        stat_gen.append(self.poolFitness())
        stat_gen.append(self.best.strandGenes)
        self.stat_board.append(stat_gen)

    def GeneticStigmergicStep(self):
        """
        One step in the Genetic Stigmergic main algorithm
        1. Updating the mating pool with current population
        2. Creating a new Generation using selection / crossover / mutation
        """
        self.updateMatingPool()
        self.newGeneration()

    def run(self):
        """
        General execution template.
        Iterates for a given number of steps
        """
        self.iteration = 0
        while self.iteration < self.maxIter:
            self.GeneticStigmergicStep()
            self.iteration += 1

        self.iprint("Total iterations: {}".format(self.iteration))
        self.iprint("Best Solution: {}".format(self.best.getStrandFitness()))

        if self.verbose == True: self.display_stats()

        self.iprint("Best Strand Details: {}".format(self.best.origGenes))
        self.iprint("Best Strand Fitness: {}".format(self.best.getStrandFitness()))
        self.iprint(self.best.strandGenes)


    def display_stats(self):

        self.iprint("[s] Display execution statistics :")
        t_iter = []
        t_slct, t_xovr, t_mutn = 0, 0, 0

        for i in range(len(self.stat_board)):
            self.iprint(
                "Gen {} - {} Sel. : {:.19f} - {} Xov. : {:.19f} - {} Mut. : {:.19f} - Best Fit : {:.9f} - Pool Fit : {:.9f}".format(
                    i + 1, self.slctType, self.stat_board[i][0], self.xovrType, self.stat_board[i][1], self.mutnType,
                    self.stat_board[i][2], self.stat_board[i][3], self.stat_board[i][4]
                ))
            t_iter.append([i + 1])
            t_slct += self.stat_board[i][0]
            t_xovr += self.stat_board[i][1]
            t_mutn += self.stat_board[i][2]

        bests_global = [stat[3] for stat in self.stat_board]
        bests_local = [stat[4] for stat in self.stat_board]

        self.iprint("[s] Timed statistics for this run :")
        self.iprint("    Total initialization time : {}".format(self.stat_inittime))
        self.iprint("    Total time for Selection : {}".format(t_slct))
        self.iprint("    Total time for Crossover : {}".format(t_xovr))
        self.iprint("    Total time for Mutations : {}".format(t_mutn))

        self.iprint("\n[v] generating plot")

        fig = mpl.figure(figsize=(12, 7), dpi=100)
        mpl.xlim(0, self.maxIter)
        mpl.plot(t_iter, bests_global)
        mpl.plot(t_iter, bests_local)
        mpl.title("Strand Fitness evol. - Pop={} - MaxI={} - Sel={} - Xov.={} - Mut.={}({})".format(self.popSize,
                  self.maxIter, self.slctType, self.xovrType, self.mutnType, self.mutationRate))
        mpl.ylabel("Strand Fitness Score")
        mpl.xlabel("Generation")
        mpl.ticklabel_format(style='plain', axis='both')

        mpl.show()
        mpl.clf()
        mpl.cla()
        mpl.close()



