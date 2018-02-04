import numpy as np
import math
import random

from deap import base
from deap import creator
from deap import tools
from deap import benchmarks
from deap.benchmarks import binary
from deap import cma
from deap import algorithms
import array
#import cma
toolbox = base.Toolbox()


# deep copy
import copy
import os.path
# delay
import time

# load and save
import pickle
import csv

# Visualization
import matplotlib.pyplot as plt

# communication
import serial # this is pyserial not serial
import struct 

# Settings or global variables
populationSize = 10
genomeLength = 15
amountTrials = 1
name = 'test'
genome = [] # template
mutationRate = 0.1
sigma = 50
useSigma = True
amountGen = 20
# not implemented yet
elitismFactor = 0.11
replaceWorst = 0.1

class Individual:
    def __init__(self):
        self.fitness = 0
        self.genome = []
 
def initializePopulation(populationSize):
    pop = []
    for i in range(populationSize):
        newGenome = []
        for j in range(genomeLength):
            newGenome.append(random.randint(0,255))  
        newInd = Individual()
        newInd.genome = newGenome
        pop.append(newInd)  
    return pop

def selectPop(pop):
    return pop

def printPop(pop):
	for i in range(len(pop)):
		print str(i) + ' has fitness ' + str(pop[i].fitness) + ' and genome ' + str(pop[i].genome) 

def selRandom(individuals, k):
    return [random.choice(individuals) for i in xrange(k)]

def selectTournament(pop, tournsize):
    k = len(pop)
    chosen = []
    for i in xrange(k):
        aspirants = selRandom(pop, tournsize)
        bestAspirant = 0
        bestAspirantFit = 0
        for j in range(len(aspirants)):
            if (aspirants[j].fitness > bestAspirantFit):
                bestAspirant = j
                bestAspirantFit = aspirants[j].fitness
        chosen.append(copy.deepcopy(aspirants[bestAspirant]))
    return chosen

def selectTournamentElitism(pop, tournsize, elitismFactor, eliteMut):
    amountElites = int(len(pop) * elitismFactor)
    elitePop = []
    eliteFitnesses = []
    eliteInds = []
    for i in range(amountElites):
        newElite = -1
        eliteFitnesses.append(0.0)
        eliteInds.append(-1)
        for j in range(len(pop)): 
            sameElite = False
            for k in range(len(eliteInds)):
                if (j == eliteInds[k]):
                    sameElite = True
            if (sameElite == False and pop[j].fitness >= eliteFitnesses[i]):
                eliteFitnesses[i] = pop[j].fitness
                eliteInds[i] = j
                newElite = j
        elitePop.append(copy.deepcopy(pop[newElite]))
    remainingPopLen = len(pop) - amountElites
    k = remainingPopLen
    chosen = []
    for i in xrange(k):
        aspirants = selRandom(pop, tournsize)
        bestAspirant = 0
        bestAspirantFit = 0
        for j in range(len(aspirants)):
            if (aspirants[j].fitness > bestAspirantFit):
                bestAspirant = j
                bestAspirantFit = aspirants[j].fitness
        chosen.append(copy.deepcopy(aspirants[bestAspirant]))
    if (eliteMut == False):
        if (useSigma == False):
            chosen = mutatePop(chosen,mutationRate) # random value
        else:
            chosen = mutatePopS(chosen,mutationRate,sigma) # with sigma distance
    elif (eliteMut == True):
        if (useSigma == False):
            chosen = mutatePop(chosen,mutationRate) # random value
            elitePop = mutatePop(elitePop,mutationRate) # random value
        else:
            chosen = mutatePopS(chosen,mutationRate,sigma) # with sigma distance
            elitePop = mutatePopS(elitePop,mutationRate,sigma) # with sigma distance
    chosenPop = elitePop + chosen
    return chosenPop


def mutatePop(pop, mr):
    for i in range(len(pop)): 
        for j in range(len(pop[i].genome)):
            if (random.uniform(0.0,1.0) < mr):
                pop[i].genome[j] = random.randint(0,255)
    return pop

def mutatePopS(pop,mr,sigma):
    for i in range(len(pop)): 
        for j in range(len(pop[i].genome)):
            if (random.uniform(0.0,1.0) < mr):
                pop[i].genome[j] += int(random.gauss(0, sigma))
                # random.randint(-sigma, +sigma)
                if (pop[i].genome[j] > 255):
                    pop[i].genome[j] = 255
                elif (pop[i].genome[j] < 0):
                    pop[i].genome[j] = 0
    return pop

def mutatePopS_B(pop,mr,sigma):
    for i in range(len(pop)): 
        for j in range(len(pop[i].genome)):
            if (random.uniform(0.0,1.0) < mr):
                pop[i].genome[j] += random.randint(-sigma, +sigma)
                if (pop[i].genome[j] > 255):
                    pop[i].genome[j] = 255
                elif (pop[i].genome[j] < 0):
                    pop[i].genome[j] = 0
    return pop


class Hof:
    def __init__(self):
        self.maxFit = []
        self.individuals = []
    def save(self, name):
        f = open(name, "w")
        pickle.dump(self.maxFit, f)
        pickle.dump(self.individuals, f)
        f.close()
    def load(self,name):
        f = open(name, "r")
        self.maxFit = pickle.load(f)
        self.individuals = pickle.load(f)
        f.close()

class Data:
    def __init__(self):
        self.individuals = []
        self.fitnessValues = []
        self.p0 = []
        self.p25 = []
        self.p50 = [] # median
        self.p75 = []
        self.p100 = []
        self.std = [] # standard deviation
        self.avg = [] # average
        self.x = []
    def save(self, name):
        f = open(name + 'data', "w")
        pickle.dump(self.fitnessValues, f)
        pickle.dump(self.p0, f)
        pickle.dump(self.p25, f)
        pickle.dump(self.p50, f)
        pickle.dump(self.p75, f)
        pickle.dump(self.p100, f)
        pickle.dump(self.std, f)
        pickle.dump(self.avg, f)
        pickle.dump(self.x, f)
        pickle.dump(self.individuals)
        f.close()
    def load(self,name):
        f = open(name + 'data', "r")
        self.fitnessValues = pickle.load(f)
        self.p0 = pickle.load(f)
        self.p25 = pickle.load(f)
        self.p50 = pickle.load(f)
        self.p75 = pickle.load(f)
        self.p100 = pickle.load(f)
        self.std = pickle.load(f)
        self.avg = pickle.load(f)
        self.x = pickle.load(f)
        self.individuals = pickle.load(f)
        f.close()
    def setValues(self,fitnessValues, individuals):
        self.fitnessValues.append(fitnessValues)
        self.p0.append(np.percentile(fitnessValues, 0))
        self.p25.append(np.percentile(fitnessValues, 25))
        self.p50.append(np.percentile(fitnessValues, 50))
        self.p75.append(np.percentile(fitnessValues, 75))
        self.p100.append(np.percentile(fitnessValues, 100))
        self.std.append(np.std(fitnessValues))
        self.avg.append(np.average(fitnessValues))
        self.x.append(len(self.x))
        self.individuals(individuals)
    def plotGraph(self, ax):
        ax.clear()
        ax.plot(self.p100, color = 'black')
        ax.plot(self.p50, color = 'white')
        ax.plot(self.p0, color = 'black')
        ax.plot(self.avg, color = 'blue')
        ax.fill_between(self.x, self.p25, self.p75, color = 'grey')

def afitFunction(genome): # just to test
    return sum(genome)

def anotherfitFunction(genome): # just to test
    return sum(genome),

def printPopulation(pop):
    for i in range(len(pop)):
        print pop[i].genome , ' has fitness: ' , pop[i].fitness

def initialPopulationEvaluation(pop, ser, data, hof):
    for i in range(populationSize):
        #print "should write : " 
        time.sleep(0.5) # pause between sending genome        
        ser.write(pop[i].genome)
        #values = bytearray(pop[i].genome)
        #values = bytearray([4, 9, 62, 144, 56, 30, 147, 3, 210, 89, 111, 78, 184, 151, 17, 129])
        #print values
        #ser.write(values)
        time.sleep(0.5) # pause between sending genome        
        print "waiting for message"
        #print ser.read()
        #print ser.readline()
        fit = float(ser.readline()) # waits for incoming line (fitness value)  
        #fit = 0.0
        print 'received robot fitness: ' , fit
        ser.flush() # remove messages in cue       
        pop[i].fitness = fit
    fitnessAr = []
    inds = []
    bestInd = 0
    bestFit = 0
    for i in range(populationSize):
        fitnessAr.append(pop[i].fitness)
        inds.append(pop[i].genome)
        if (pop[i].fitness > bestFit):
            bestFit = pop[i].fitness
            bestInd = i
    data.setValues(fitnessAr, inds)
    hof.individuals.append(pop[bestInd])
    hof.maxFit.append(bestFit)
    printPop(pop)

def initialPopulationEvaluationTest(pop, data, hof):
    for i in range(populationSize):
        fit = afitFunction(pop[i].genome)
        pop[i].fitness = fit
        print "fitness: " , fit 
    fitnessAr = []
    bestInd = 0
    bestFit = 0
    inds = []
    for i in range(populationSize):
        fitnessAr.append(pop[i].fitness)
        inds.append(pop[i].genome)
        if (pop[i].fitness > bestFit):
            bestFit = pop[i].fitness
            bestInd = i
    data.setValues(fitnessAr, inds)
    hof.individuals.append(pop[bestInd])
    hof.maxFit.append(bestFit)
    return pop

def continueEvolution(pop, hof, ax, data, ser):
    ser.flush()
    time.sleep(0.5)
	# elitism tournament includes mutation
    pop = selectTournamentElitism(pop,3, elitismFactor, False)
    for i in range(populationSize):
        fits_t = [] # stores fitness values of all individuals
        for t in range(amountTrials):
            time.sleep(0.5)
            print "sending genome: " + str(pop[i].genome)
            ser.write(pop[i].genome)
            time.sleep(0.5)
            fit = float(ser.readline()) # waits until line is received
            print 'received robot fitness: ', fit
            fits_t.append(fit)
            ser.flush()
        listSum = sum(fits_t)
        listLength = len(fits_t)
        listAverage = listSum / listLength
        pop[i].fitness = listAverage # lowest value will be fitness 
        print "pop fitness of ind "+ str(i) + " = " + str(pop[i].fitness)    
    bestInd = 0
    bestFit = 0
    fitnessAr = []
    inds = []
    for i in range(populationSize):
        fitnessAr.append(pop[i].fitness)
        inds.append(pop[i].genome)
        if (pop[i].fitness > bestFit):
            bestFit = pop[i].fitness
            bestInd = i
        
    data.setValues(fitnessAr, inds)
    hof.individuals.append(pop[bestInd])
    hof.maxFit.append(bestFit)
#    ax.clear()
#    ax.plot(hof.maxFit)
    data.plotGraph(ax)
    printPop(pop)
    plt.pause(0.05)
    return pop


def continueEvolutionTest(pop, hof, ax, data):
    time.sleep(0.001)
#    pop = selectTournament(pop, 3)
	# elitism tournament includes mutation
    pop = selectTournamentElitism(pop,3, elitismFactor, False)
    for i in range(populationSize):
        fits_t = [] # stores fitness values of all individuals
        for t in range(amountTrials):
            time.sleep(0.001)
            print "testing genome: " + str(pop[i].genome)
            time.sleep(0.001)
            fit = afitFunction(pop[i].genome) # waits until line is received
            print 'fitness of genome ', i , ' : ', fit
            fits_t.append(fit)
        pop[i].fitness = min(fits_t) # lowest value will be fitness 
        print "pop fitness of ind "+ str(i) +" = " + str(pop[i].fitness)    
    bestInd = 0
    bestFit = 0
    fitnessAr = []
    inds = []
    for i in range(populationSize):
        fitnessAr.append(pop[i].fitness)
        inds.append(pop[i].genome)
        if (pop[i].fitness > bestFit):
             bestFit = pop[i].fitness
             bestInd = i
     
    data.setValues(fitnessAr, inds)
    hof.individuals.append(pop[bestInd])
    hof.maxFit.append(bestFit)
#    ax.clear()
    data.plotGraph(ax)
#    ax.plot(hof.maxFit, color = 'black')
    plt.pause(0.001)
    return pop

def evaluateIndividual(genome, ser):
    time.sleep(0.5)
    print "sending genome: " + str(genome)
    ser.write(pop[i].genome)
    time.sleep(0.5)
    fit = float(ser.readline()) # waits until line is received
    print 'received robot fitness: ', fit
    fits_t.append(fit)
    ser.flush()
    return fitness

def run(hof, data, ax, name):
    maxFit = 0
    # open serial port
    ser = serial.Serial('COM21', 9600)
    print "opened serial"
    time.sleep(0.5) # pause between sending genome        
    ser.read_all() # necessary for getting rid of the signals used to establish communication

    ax1 = ax
    
    # create initial population of random individuals 
    pop = initializePopulation(populationSize)

    # test arduino
    # genomeTest = [120, 120, 120]
    # time.sleep(0.5)
    # plt.pause(0.05)

    # evaluate initial population
    #print "evaluating initial population"
    #initialPopulationEvaluation(pop, ser, data, hof)
    #hof.save(name=name)
    #data.save(name=name)

    # for testing without serial:
    #initialPopulationEvaluationTest(pop, data, hof)

    # CMA
    process = 1
    if (process == 1):
        N = genomeLength
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMax)
        creator.create("Strategy", array.array, typecode="d")

        def initES(icls, scls, size, imin, imax, smin, smax):
            ind = icls(random.uniform(imin, imax) for _ in range(size))
            ind.strategy = scls(random.uniform(smin, smax) for _ in range(size))
            return ind

        IND_SIZE=genomeLength
    
        MIN_VALUE, MAX_VALUE = -500., 500.
        MIN_STRAT, MAX_STRAT = -1., 1. 
        AGE = 0

        toolbox.register("attr_float", random.random)
        toolbox.register("individual", initES, creator.Individual,
                     creator.Strategy, IND_SIZE, MIN_VALUE, MAX_VALUE, MIN_STRAT, 
                     MAX_STRAT)

        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=sigma, indpb=0.2)
        toolbox.register("select", tools.selTournament, tournsize=3)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        finalPops = []
        allAllMax = []
        allAllMin = []
        finalPops = []
        finalDistr = []
        finalFits = []
        finalBestPops = []

        strategy = cma.Strategy(centroid=[125]*N, sigma=sigma, lambda_=populationSize)# lambda_=20*N)
        toolbox.register("generate",strategy.generate,creator.Individual)
        toolbox.register("update",strategy.update)
        pop = toolbox.population(n=populationSize) 
        sigmaValues = []
        fPop = []
        allMin = []
        fits = []

        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)
        logbook = tools.Logbook()
        logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])
        
        for i in range(amountGen):
            fitnesses = []
            population = toolbox.generate()
            for k in range(len(population)):
                for l in range(genomeLength):
                    if (population[k][l] > 255):
                        population[k][l] = 255
                    if (population[k][l] < 0):
                        population[k][l] = 0
                # Evaluate the individuals
                # test:
                # fitnesses.append(anotherfitFunction(population[k]))
                # with serial:
                fitnesses.append(evaluateIndividual(population[k],ser))
            for ind, fit in zip(population, fitnesses):
                ind.fitness.values = fit
            minFit = min(fitnesses)
            allMin.append(minFit[0])
            # Update the strategy with the evaluated individuals
            toolbox.update(population)
            sigmaValues.append(strategy.sigma)
            # Update the hall of fame and the statistics with the
            # currently evaluated population
            # hof.update(population)
            record = stats.compile(population)
            logbook.record(evals=len(population), gen=i, **record)
            print logbook.stream
            fPop.append(population)

            bestInd = 0
            bestFit = 0
            fitnessAr = []
            inds = []
            for i in range(populationSize):
                fitnessAr.append(fitnesses[i])
                inds.append(population[i])
                if (fitnesses[i] > bestFit):
                    bestFit = fitnesses[i]
                    bestInd = i
        
            data.setValues(fitnessAr, inds)
            theInd = Individual()
            theInd.genome = population[bestInd]
            theInd.fitness = fitnesses[i]
            hof.individuals.append(theInd)
            hof.maxFit.append(bestFit)
#    ax.clear()
#    ax.plot(hof.maxFit)
        data.plotGraph(ax)
        finalPops.append(fPop)
        allAllMin.append(allMin)
        ax1.plot(allMin)
       
       
        #pop = initializePopulation(populationSize)
        #pop = initialPopulationEvaluationTest(pop,data, hof)
        #strategy = cma.Strategy(centroid=[125]*N, sigma=sigma, lambda_=populationSize)
        #for i in range(10):
        #    pop = strategy.update(pop, hof, ax1, data)
        #    pop = initialPopulationEvaluationTest(pop,data, hof)
        
    # evolutionary loop
    print "Starting Evolution"
    #for i in range(amountGen):
    #    pop = continueEvolution(pop, hof, ax1, data, ser)
    #    hof.save(name=name)
    #    data.save(name=name)
        

    #for i in range(amountGen):
    #    pop = continueEvolutionTest(pop, hof, ax1, data)
    data.plotGraph(ax1) 
    #plt.show()
    return hof, data


def main():
	
    amountRuns = 2
    hofs = []
    datas = []
    plt.ion()
    fig = plt.figure()
    name = 'ExpZero'
    for i in range(amountRuns):
        counter = 0
        nameFound = False 
        while (nameFound == False):
            fname = name + str(counter) + '_' + '0'
            if (os.path.isfile(fname) == True):
                counter+=1
                fname = name + str(counter) + '_' + '0'
                nameFound = False
            else:
                name = name + str(counter) + '_'
                nameFound = True
        x = amountRuns
        ax = fig.add_subplot(x,1,i + 1)    
        hof = Hof()
        data = Data()
        run(hof,data, ax, name + str(i))
        hofs.append(hof)
        datas.append(data)
    plt.ioff()
    plotDatas(datas)
    plt.show()


def plotDatas(datas):
    fig = plt.figure()
    max100 = []
    max0 = []
    max50 = []
    max25 = []
    max75 = []
    x = []
    for n in range(len(datas[0].p100)):
        fits = []
        for i in range(len(datas)):
            fits.append(datas[i].p100[n])
        max100.append(np.percentile(fits,100))
        max0.append(np.percentile(fits,0))
        max50.append(np.percentile(fits,50))
        max25.append(np.percentile(fits,25))
        max75.append(np.percentile(fits,75))
        x.append(len(x))
    ax = fig.add_subplot(1,1,1)
    ax.plot(max100, color = 'black') 
    ax.plot(max0, color = 'black') 
    ax.plot(max50, color = 'white') 
    ax.fill_between(x,max25,max75, color = 'grey')

def loadBestIndividual():
    filename = 'test0_0'
    hof=Hof()
    hof.load(filename)
    bestIndividual = []
    bestGen = 0
    mf = 0
    ser = serial.Serial('COM21', 9600)
    time.sleep(0.5) # pause between sending genome        
    ser.read_all()
    time.sleep(0.5) # pause between sending genome        
    
    for i in range(len(hof.maxFit)):
        if (hof.maxFit[i] > mf):
            bestIndividual = hof.individuals[i].genome
            mf = hof.maxFit[i]
            bestGen = i
    print "best individual was from generation " , bestGen, " and had genome " , bestIndividual , " with a fitness of " , mf
    time.sleep(1.0) # pause between sending genome        
    ser.read_all()
    ser.write(bestIndividual)
    time.sleep(0.5) # pause between sending genome        
    print "waiting for message"
    fit = float(ser.readline()) # waits for incoming line (fitness value)  
    print 'received robot fitness: ' , fit

def loadBestIndividual(fname):
    filename = fname
    hof=Hof()
    hof.load(filename)
    bestIndividual = []
    bestGen = 0
    mf = 0
    ser = serial.Serial('COM21', 9600)
    time.sleep(0.5) # pause between sending genome        
    ser.read_all()
    time.sleep(0.5) # pause between sending genome        
    
    for i in range(len(hof.maxFit)):
        if (hof.maxFit[i] > mf):
            bestIndividual = hof.individuals[i].genome
            mf = hof.maxFit[i]
            bestGen = i
    print "best individual was from generation " , bestGen, " and had genome " , bestIndividual , " with a fitness of " , mf
    time.sleep(1.0) # pause between sending genome        
    ser.read_all()
    ser.write(bestIndividual)
    time.sleep(0.5) # pause between sending genome        
    print "waiting for message"
    fit = float(ser.readline()) # waits for incoming line (fitness value)  
    print 'received robot fitness: ' , fit


def getServoPL(bestIndividual): 
    servoPlots = []
    maxServoRange = 40
    steps = 100
    for i in range(6):
        servoPL = []
        f = 0
        for x in range(steps):
            sum = 0 
            for n in range(0,len(bestIndividual),3): 
                AMP = bestIndividual[n] / 255.0 * maxServoRange
                PHASE = bestIndividual[n + 1] / 255.0 * 4
                FREQ = bestIndividual[n + 2] / 255.0 / 5 
                f = AMP * math.sin((math.pi * ((-x) * FREQ)) + (PHASE * i))
                if (len(bestIndividual) <= 3):
                    sum = sum + f ;
                else:
                    sum = sum + (f * 2);
                if (sum > maxServoRange): sum = maxServoRange
                if (sum < -maxServoRange): sum = -maxServoRange
            servoPL.append(sum)
        servoPlots.append(servoPL)
    return servoPlots


def loadAll(fname):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)    
    hof=Hof()
    data = Data()
    hof.load(fname)
    data.load(fname)
    data.plotGraph(ax)
    hofAMP = [] 
    hofPHASE = []
    hofFREQ = []
    if (len(hof.individuals[0].genome) > 3):
        for i in range(len(hof.individuals)):
            amp = []
            phase = []
            freq = []
            for n in range(5):
                amp.append(hof.individuals[i].genome[(n * 3) + 0]) 
                phase.append(hof.individuals[i].genome[(n * 3) + 1]) 
                freq.append(hof.individuals[i].genome[(n * 3) + 2]) 
            hofAMP.append(amp) 
            hofPHASE.append(phase) 
            hofFREQ.append(freq) 

        ax2.plot(hofAMP, color = 'black')
        ax2.plot(hofPHASE, color = 'green')
        ax2.plot(hofFREQ, color = 'blue')
    else: 
        for i in range(len(hof.individuals)):
            hofAMP.append(hof.individuals[i].genome[0]) 
            hofPHASE.append(hof.individuals[i].genome[1]) 
            hofFREQ.append(hof.individuals[i].genome[2]) 
        ax2.plot(hofAMP)
        ax2.plot(hofPHASE)
        ax2.plot(hofFREQ)
    # plot HOF best function
    bestIndividual = []
    bestGen = 0
    mf = 0
    for i in range(len(hof.maxFit)):
        if (hof.maxFit[i] > mf):
            bestIndividual = hof.individuals[i].genome
            mf = hof.maxFit[i]
            bestGen = i
    print "best individual was from generation " , bestGen, " and had genome " , bestIndividual , " with a fitness of " , mf
    
    fig3 = 0
    axarr =  0 
    fig5 = plt.figure()
    showTwo = True
    showFourier = True
    ax3 = 0
    a3man = 0
    if (showTwo == True and showFourier != True):
        fig3, axarr = plt.subplots(3, sharex=True)
        axarr[2].set_xlabel('Time in milliseconds')

        #ax3 = fig3.add_subplot(3,1,2)  
        #ax3man = fig3.add_subplot(3,1,1)  
    elif (showFourier == True):
        fig3, axarr = plt.subplots(2, sharex=True)
        axarr[1].set_xlabel('Time in milliseconds')
    else:
        fig3, axarr = plt.subplots(2, sharex=True)
        axarr[1].set_xlabel('Time in milliseconds')
    
    if (showFourier == False):
        axarr[1].set_title("Run 1")
        axarr[0].set_title("Manual")
        xVals = np.arange(0, 100 * 20, 20)
        bestIndividualMan = [255,64,100];
        servoPlotsMan = getServoPL(bestIndividualMan)
        axarr[0].plot(xVals, servoPlotsMan[0])
        axarr[0].plot(xVals,servoPlotsMan[1], '--')
        axarr[0].plot(xVals,servoPlotsMan[2], '--')
        servoPlots = getServoPL(bestIndividual)
        axarr[1].plot(xVals, servoPlots[0])
        axarr[1].plot(xVals,servoPlots[1], '--')
        axarr[1].plot(xVals,servoPlots[2], '--')
        fig3.text(0.04, 0.5, 'Servo Angle', va='center', rotation='vertical')
    else:
        axarr[1].set_title("Run 2")
        axarr[0].set_title("Run 1")
        xVals = np.arange(0, 100 * 20, 20)
        servoPlots = getServoPL(bestIndividual)
        axarr[0].plot(xVals, servoPlots[0])
        axarr[0].plot(xVals,servoPlots[1], '--')
        axarr[0].plot(xVals,servoPlots[2], '--')
        fig3.text(0.04, 0.5, 'Servo Angle', va='center', rotation='vertical')

        

    ax5 = fig5.add_subplot(1,1,1)
    ax5.plot(hof.maxFit, label = 'Run 1')
      
    

    ax5.set_xlabel('Generation')
    ax5.set_ylabel('Fitness Value')
    
    if (showTwo == True): 
        
        hof2=Hof()
        hof2.load(fname + "_1")
        ax5.plot(hof2.maxFit, label = 'Run 2')
        bestIndividual2 = []
        bestGen2 = 0
        mf2 = 0
        for i in range(len(hof2.maxFit)):
            if (hof2.maxFit[i] > mf2):
                bestIndividual2 = hof2.individuals[i].genome
                mf2 = hof2.maxFit[i]
                bestGen2 = i
        print "best individual was from generation " , bestGen2, " and had genome " , bestIndividual2, " with a fitness of " , mf2
        servoPlots2 = getServoPL(bestIndividual2)
        if (showFourier == False):
            axarr[2].set_title("Run 2")
            axarr[2].plot(xVals, servoPlots2[0])
            axarr[2].plot(xVals,servoPlots2[1], '--')
            axarr[2].plot(xVals,servoPlots2[2], '--')
        else:
            axarr[1].set_title("Run 2")
            axarr[1].plot(xVals, servoPlots2[0])
            axarr[1].plot(xVals,servoPlots2[1], '--')
            axarr[1].plot(xVals,servoPlots2[2], '--')
    
    handles, labels = ax5.get_legend_handles_labels()
    
    ax5.legend(handles, labels,loc='lower right')
#    ax4.plot(xVals, servoPlots2[3], '--')
#    ax4.plot(xVals, servoPlots2[4], '--')
#    ax4.plot(xVals, servoPlots2[5], '--')
      #sum = 90.0 + (sum + 0.5); //0.5 added to round off correctly when converting to byte below
      #wave[j][i] = (byte)sum;

#fname = 'ExpZero42_0' 
#loadAll(fname)
main()
#loadBestIndividual(fname)
#plt.show()
# main()
#loadBestIndividual()
