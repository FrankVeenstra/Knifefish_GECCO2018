import numpy as np
import math
import random

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
genomeLength = 3
amountTrials = 1
name = 'test'
genome = [] # template
mutationRate = 0.2
sigma = 100
useSigma = True
amountGen = 100
# not implemented yet
elitism = 0.1
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

def mutatePop(pop, mr):
    for i in range(len(pop)): 
        for j in range(len(pop[i].genome)):
            if (random.uniform(0,1) < mr):
                pop[i].genome[j] = random.randint(0,255)
    return pop

def mutatePopS(pop,mr,sigma):
    for i in range(len(pop)): 
        for j in range(len(pop[i].genome)):
            if (random.uniform(0,1) < mr):
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
        f.close()
    def setValues(self,fitnessValues):
        self.fitnessValues.append(fitnessValues)
        self.p0.append(np.percentile(fitnessValues, 0))
        self.p25.append(np.percentile(fitnessValues, 25))
        self.p50.append(np.percentile(fitnessValues, 50))
        self.p75.append(np.percentile(fitnessValues, 75))
        self.p100.append(np.percentile(fitnessValues, 100))
        self.std.append(np.std(fitnessValues))
        self.avg.append(np.average(fitnessValues))
        self.x.append(len(self.x))
    def plotGraph(self, ax):
        ax.clear()
        ax.plot(self.p100, color = 'black')
        ax.plot(self.p50, color = 'white')
        ax.plot(self.p0, color = 'black')
        ax.plot(self.avg, color = 'blue')
        ax.fill_between(self.x, self.p25, self.p75, color = 'grey')

def afitFunction(genome): # just to test
    return sum(genome)

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
    bestInd = 0
    bestFit = 0
    for i in range(populationSize):
        fitnessAr.append(pop[i].fitness)
        if (pop[i].fitness > bestFit):
            bestFit = pop[i].fitness
            bestInd = i
    data.setValues(fitnessAr)
    hof.individuals.append(pop[bestInd])
    hof.maxFit.append(bestFit)
    print(pop)

def initialPopulationEvaluationTest(pop, data, hof):
    for i in range(populationSize):
        fit = afitFunction(pop[i].genome)
        pop[i].fitness = fit
        print "fitness: " , fit 
    fitnessAr = []
    bestInd = 0
    bestFit = 0
    for i in range(populationSize):
        fitnessAr.append(pop[i].fitness)
        if (pop[i].fitness > bestFit):
            bestFit = pop[i].fitness
            bestInd = i
    data.setValues(fitnessAr)
    hof.individuals.append(pop[bestInd])
    hof.maxFit.append(bestFit)

def continueEvolution(pop, hof, ax, data, ser):
    ser.flush()
    time.sleep(0.5)
    pop = selectTournament(pop, 3)
    if (useSigma == False):
        pop = mutatePop(pop,mutationRate) # random value
    else:
        pop = mutatePopS(pop,mutationRate,sigma) # with sigma distance
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
        pop[i].fitness = min(fits_t) # lowest value will be fitness 
        print "pop fitness of ind "+ str(i) + " = " + str(pop[i].fitness)    
    bestInd = 0
    bestFit = 0
    fitnessAr = []

    for i in range(populationSize):
        fitnessAr.append(pop[i].fitness)
        if (pop[i].fitness > bestFit):
            bestFit = pop[i].fitness
            bestInd = i
        
    data.setValues(fitnessAr)
    hof.individuals.append(pop[bestInd])
    hof.maxFit.append(bestFit)
#    ax.clear()
#    ax.plot(hof.maxFit)
    data.plotGraph(ax)
    sum = 0
    print(pop)
    plt.pause(0.05)


def continueEvolutionTest(pop, hof, ax, data):
    time.sleep(0.001)
    pop = selectTournament(pop, 3)
    if (useSigma == False):
        pop = mutatePop(pop,mutationRate) # random value
    else:
        pop = mutatePopS(pop,mutationRate,sigma) # with sigma distance
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
    for i in range(populationSize):
        fitnessAr.append(pop[i].fitness)
        if (pop[i].fitness > bestFit):
             bestFit = pop[i].fitness
             bestInd = i
     
    data.setValues(fitnessAr)
    hof.individuals.append(pop[bestInd])
    hof.maxFit.append(bestFit)
#    ax.clear()
    data.plotGraph(ax)
#    ax.plot(hof.maxFit, color = 'black')
    sum = 0
    plt.pause(0.001)
    return pop

def run(hof, data, ax, name):
    maxFit = 0
    # open serial port
    ser = serial.Serial('COM21', 9600)
    print "opened serial"
    time.sleep(0.5) # pause between sending genome        
    ser.read_all()

    ax1 = ax
    
    # create initial population of random individuals 
    pop = initializePopulation(populationSize)

    # test arduino
    # genomeTest = [120, 120, 120]
    # time.sleep(0.5)
    # plt.pause(0.05)

    # evaluate initial population
    print "evaluating initial population"
    initialPopulationEvaluation(pop, ser, data, hof)
    hof.save(name=name)
    data.save(name=name)
        

    # for testing without serial:
    #initialPopulationEvaluationTest(pop, data, hof)
    
    # evolutionary loop
    print "Starting Evolution"
    for i in range(amountGen):
        continueEvolution(pop, hof, ax1, data, ser)
        hof.save(name=name)
        data.save(name=name)
        

    #for i in range(amountGen):
    #    pop = continueEvolutionTest(pop, hof, ax1, data)
    data.plotGraph(ax1) 
    #plt.show()
    return hof, data


def main():
	
    amountRuns = 1
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

def loadAll(fname):
    fig = plt.figure()
    ax = fig.add_subplot(1,1, 1)    
    hof=Hof()
    data = Data()
    hof.load(fname)
    data.load(fname)
    data.plotGraph(ax) 

fname = 'ExpZero6_0' 
loadAll(fname)
#main()
loadBestIndividual(fname)
plt.show()
