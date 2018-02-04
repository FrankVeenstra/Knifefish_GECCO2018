#    This file is part of DEAP.
#
#    DEAP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    DEAP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with DEAP. If not, see <http://www.gnu.org/licenses/>.

#    Special thanks to Nikolaus Hansen for providing major part of
#    this code. The CMA-ES algorithm is provided in many other languages
#    and advanced versions at http://www.lri.fr/~hansen/cmaesintro.html.

"""A module that provides support for the Covariance Matrix Adaptation
Evolution Strategy.
"""
import copy
from math import sqrt, log, exp
import numpy

class Strategy(object):
    """
    A strategy that will keep track of the basic parameters of the CMA-ES
    algorithm ([Hansen2001]_).

    :param centroid: An iterable object that indicates where to start the
                     evolution.
    :param sigma: The initial standard deviation of the distribution.
    :param parameter: One or more parameter to pass to the strategy as
                      described in the following table, optional.

    +----------------+---------------------------+----------------------------+
    | Parameter      | Default                   | Details                    |
    +================+===========================+============================+
    | ``lambda_``    | ``int(4 + 3 * log(N))``   | Number of children to      |
    |                |                           | produce at each generation,|
    |                |                           | ``N`` is the individual's  |
    |                |                           | size (integer).            |
    +----------------+---------------------------+----------------------------+
    | ``mu``         | ``int(lambda_ / 2)``      | The number of parents to   |
    |                |                           | keep from the              |
    |                |                           | lambda children (integer). |
    +----------------+---------------------------+----------------------------+
    | ``cmatrix``    | ``identity(N)``           | The initial covariance     |
    |                |                           | matrix of the distribution |
    |                |                           | that will be sampled.      |
    +----------------+---------------------------+----------------------------+
    | ``weights``    | ``"superlinear"``         | Decrease speed, can be     |
    |                |                           | ``"superlinear"``,         |
    |                |                           | ``"linear"`` or            |
    |                |                           | ``"equal"``.               |
    +----------------+---------------------------+----------------------------+
    | ``cs``         | ``(mueff + 2) /           | Cumulation constant for    |
    |                | (N + mueff + 3)``         | step-size.                 |
    +----------------+---------------------------+----------------------------+
    | ``damps``      | ``1 + 2 * max(0, sqrt((   | Damping for step-size.     |
    |                | mueff - 1) / (N + 1)) - 1)|                            |
    |                | + cs``                    |                            |
    +----------------+---------------------------+----------------------------+
    | ``ccum``       | ``4 / (N + 4)``           | Cumulation constant for    |
    |                |                           | covariance matrix.         |
    +----------------+---------------------------+----------------------------+
    | ``ccov1``      | ``2 / ((N + 1.3)^2 +      | Learning rate for rank-one |
    |                | mueff)``                  | update.                    |
    +----------------+---------------------------+----------------------------+
    | ``ccovmu``     | ``2 * (mueff - 2 + 1 /    | Learning rate for rank-mu  |
    |                | mueff) / ((N + 2)^2 +     | update.                    |
    |                | mueff)``                  |                            |
    +----------------+---------------------------+----------------------------+

    .. [Hansen2001] Hansen and Ostermeier, 2001. Completely Derandomized
       Self-Adaptation in Evolution Strategies. *Evolutionary Computation*

    """
    def __init__(self, centroid, sigma, **kargs):
        self.params = kargs

        # Create a centroid as a numpy array
        self.centroid = numpy.array(centroid)

        self.dim = len(self.centroid)
        self.sigma = sigma
        self.pc = numpy.zeros(self.dim)
        self.ps = numpy.zeros(self.dim)
        self.chiN = sqrt(self.dim) * (1 - 1. / (4. * self.dim) +
                                      1. / (21. * self.dim ** 2))

        self.C = self.params.get("cmatrix", numpy.identity(self.dim))
        self.diagD, self.B = numpy.linalg.eigh(self.C)

        indx = numpy.argsort(self.diagD)
        self.diagD = self.diagD[indx] ** 0.5
        self.B = self.B[:, indx]
        self.BD = self.B * self.diagD

        self.cond = self.diagD[indx[-1]] / self.diagD[indx[0]]

        self.lambda_ = self.params.get("lambda_", int(4 + 3 * log(self.dim)))
        self.update_count = 0
        self.computeParams(self.params)

    def generate(self, ind_init):
        """Generate a population of :math:`\lambda` individuals of type
        *ind_init* from the current strategy.

        :param ind_init: A function object that is able to initialize an
                         individual from a list.
        :returns: A list of individuals.
        """
        arz = numpy.random.standard_normal((self.lambda_, self.dim))
        arz = self.centroid + self.sigma * numpy.dot(arz, self.BD.T)
        return map(ind_init, arz)

    def update(self, population, hof, ax, data):
        """Update the current covariance matrix strategy from the
        *population*.

        :param population: A list of individuals from which to update the
                           parameters.
        """
        population.sort(key=lambda ind: ind.fitness, reverse=True)

        popGenomes = []
        for i in range(len(population)):
            genome = []
            for j in range(len(population[i].genome)):
                if (population[i].genome[j] > 255):
                    population[i].genome[j] = 255
                elif (population[i].genome[j] < 0):
                    population[i].genome[j] = 0
                genome.append(population[i].genome[j])
            popGenomes.append(genome)

        old_centroid = self.centroid
        self.centroid = numpy.dot(self.weights, popGenomes[0:self.mu])

        c_diff = self.centroid - old_centroid

        # Cumulation : update evolution path
        self.ps = (1 - self.cs) * self.ps \
            + sqrt(self.cs * (2 - self.cs) * self.mueff) / self.sigma \
            * numpy.dot(self.B, (1. / self.diagD) *
                        numpy.dot(self.B.T, c_diff))

        hsig = float((numpy.linalg.norm(self.ps) /
                      sqrt(1. - (1. - self.cs) ** (2. * (self.update_count + 1.))) / self.chiN <
                      (1.4 + 2. / (self.dim + 1.))))

        self.update_count += 1

        self.pc = (1 - self.cc) * self.pc + hsig \
            * sqrt(self.cc * (2 - self.cc) * self.mueff) / self.sigma \
            * c_diff

        # Update covariance matrix
        artmp = popGenomes[0:self.mu] - old_centroid
        self.C = (1 - self.ccov1 - self.ccovmu + (1 - hsig) *
                  self.ccov1 * self.cc * (2 - self.cc)) * self.C \
            + self.ccov1 * numpy.outer(self.pc, self.pc) \
            + self.ccovmu * numpy.dot((self.weights * artmp.T), artmp) \
            / self.sigma ** 2

        self.sigma *= numpy.exp((numpy.linalg.norm(self.ps) / self.chiN - 1.) *
                                self.cs / self.damps)

        self.diagD, self.B = numpy.linalg.eigh(self.C)
        indx = numpy.argsort(self.diagD)

        self.cond = self.diagD[indx[-1]] / self.diagD[indx[0]]

        self.diagD = self.diagD[indx] ** 0.5
        self.B = self.B[:, indx]
        self.BD = self.B * self.diagD

        for i in range(len(population)):
            for j in range(len(population[i].genome)):
                if (popGenomes[i][j] > 255):
                    popGenomes[i][j] = 255
                elif (popGenomes[i][j] < 0):
                    popGenomes[i][j] = 0
                population[i].genome[j] = popGenomes[i][j]
            population[i].fitness = 0
        return population

    def computeParams(self, params):
        """Computes the parameters depending on :math:`\lambda`. It needs to
        be called again if :math:`\lambda` changes during evolution.

        :param params: A dictionary of the manually set parameters.
        """
        self.mu = params.get("mu", int(self.lambda_ / 2))
        rweights = params.get("weights", "superlinear")
        if rweights == "superlinear":
            self.weights = log(self.mu + 0.5) - \
                numpy.log(numpy.arange(1, self.mu + 1))
        elif rweights == "linear":
            self.weights = self.mu + 0.5 - numpy.arange(1, self.mu + 1)
        elif rweights == "equal":
            self.weights = numpy.ones(self.mu)
        else:
            raise RuntimeError("Unknown weights : %s" % rweights)

        self.weights /= sum(self.weights)
        self.mueff = 1. / sum(self.weights ** 2)

        self.cc = params.get("ccum", 4. / (self.dim + 4.))
        self.cs = params.get("cs", (self.mueff + 2.) /
                             (self.dim + self.mueff + 3.))
        self.ccov1 = params.get("ccov1", 2. / ((self.dim + 1.3) ** 2 +
                                               self.mueff))
        self.ccovmu = params.get("ccovmu", 2. * (self.mueff - 2. +
                                                 1. / self.mueff) /
                                 ((self.dim + 2.) ** 2 + self.mueff))
        self.ccovmu = min(1 - self.ccov1, self.ccovmu)
        self.damps = 1. + 2. * max(0, sqrt((self.mueff - 1.) /
                                           (self.dim + 1.)) - 1.) + self.cs
        self.damps = params.get("damps", self.damps)
