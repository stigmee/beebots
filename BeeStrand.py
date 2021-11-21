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
## File : BeeStrand.py
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
#   Add a stigmergic evaluation function

import random

class Strand:

    def __init__(self, _strandSize, _origGenes, _childGenes, _nni=False):
        """
        Genetic setup of a strand
        To be used during testing. Fitness evaluation is based on the best ascending order
        compared to the original indexing. A way to measure the fitness is by inversion measurement
        https://en.wikipedia.org/wiki/Inversion_(discrete_mathematics)
        """
        self.strandFitness = 0
        self.strandGenes = []
        self.strandSize = _strandSize
        # Orig genes is a given dictionary of (key, values) tuples which can be used for initialization
        self.origGenes = _origGenes

        if _childGenes:
            # Genes resulting from a crossover operation, not to be shuffled
            self.strandGenes = _childGenes
        else:
            if _nni == False:
                # If we want random initialization, we just shuffle the incoming data
                self.strandGenes = list(_origGenes.keys())
                random.shuffle(self.strandGenes)
            else:
                # If we want NNI we keep the incoming indexing and just chose a random starting point
                i = random.randint(0,self.strandSize-1)
                self.strandGenes[i:self.strandSize-1] = _origGenes.keys()[i:_origGenes-1]
                self.strandGenes[0:i-1] = _origGenes.keys()[0:i-1]
        # Compute the fitness of a strand
        self.computeStrandFitness()

    def copy(self):
        """
        duplicating a strand (in case we intend to copy the same one for partial population replacement)
        """
        s = Strand(self.strandSize, self.origGenes, self.strandGenes[0:self.strandSize])
        s.strandFitness = self.getStrandFitness()
        return s

    def computeStrandFitness(self):
        """
        The inversion fitness of a strand computes how well this strand is sorted,
        given the original sort index and the actual strand genes (url indexing of the strand instance)
        We compute the number of inversions for all possible pairs of indices and add the initial index value

        For testing purpose, if we want to evaluate a strand originally indexes,
        we try to minimize this value. This computation does not apply to user curation of strands
        """
        distances = 0
        for i in range(0, self.strandSize - 1):
            # Add the distance to 0 for each gene
            distances += self.strandGenes[i]
            for j in range(i+1, self.strandSize - 1):
                if self.strandGenes[i] > self.strandGenes[j]:
                    # If a permutation is detected, we add the corresponding index distance
                    distances += self.strandGenes[i] - self.strandGenes[j]
        self.strandFitness = distances


    def getStrandFitness(self):
        """
        return the fitness of a strand.
        """
        return self.strandFitness
