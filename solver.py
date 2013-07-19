#!/usr/bin/python
# -*- coding: utf-8 -*-
import datetime


def trivialSolution(nodeCount):
    # build a trivial solution
    # every node has its own color
    solution = range(0, nodeCount)

    # prepare the solution in the specified output format
    outputData = str(nodeCount) + ' ' + str(0) + '\n'
    outputData += ' '.join(map(str, solution))

    return outputData

def secondTrivial(nodeCount, edges):
    #Create tuples of form (node, degree)
    degree_tuples = []
    for node in range(0,nodeCount):
        nodeEdges = [[x,y] for [x,y] in edges if x==node or y==node]
        degree = len(nodeEdges)
        degree_tuples.append([node, degree])

    #Sort by Degree
    degree_tuples=sorted(degree_tuples, key=lambda node: node[1], reverse=True)

    #Create tuples of form (node, degree, color)
    color_tuples = []
    i=0
    for degree_tuple in degree_tuples:
        color_tuples.append([degree_tuple[0],degree_tuple[1],-1])
        i=i+1

    #use heuristic attempt algorithm to solve
    currentColor = 0
    #color_tuples[0][2]=0
    while len([x for [x,y,z] in color_tuples if z ==-1])>0:
        initial_node_index = color_tuples.index([[x,y,z] for [x,y,z] in color_tuples if z ==-1][0])
        color_tuples[initial_node_index][2]=currentColor
        if len([x for [x,y,z] in color_tuples if z ==-1])>0:
            for node in color_tuples:
                node_index = color_tuples.index([[x,y,z] for [x,y,z] in color_tuples if x ==node[0]][0])
                if(node[2]==-1):
                    nodeEdges = [[x,y] for [x,y] in edges if x==node[0] or y==node[0]]
                    connections = {a for [a,b] in nodeEdges if b==node[0]}|{d for [c,d] in nodeEdges if c==node[0]}
                    currentColorNodes = {x for [x,y,z] in color_tuples if z==currentColor}
                    if len(connections&currentColorNodes)== 0:
                        color_tuples[node_index][2]=currentColor
            currentColor += 1

    # prepare the solution in the specified output format
    numColors = max([z for (x,y,z) in color_tuples])+1
    solution = [z for (x,y,z) in sorted(color_tuples, key=lambda node: node[0])]
    outputData = str(numColors) + ' ' + str(0) + '\n'
    outputData += ' '.join(map(str, solution))

    return outputData

def bruteTrivial(nodeCount, edges, offset):
    #Create tuples of form (node, degree)
    degree_tuples = []
    for node in range(0,nodeCount):
        nodeEdges = [[x,y] for [x,y] in edges if x==node or y==node]
        degree = len(nodeEdges)
        degree_tuples.append([node, degree])

    #Sort by Degree
    degree_tuples=sorted(degree_tuples, key=lambda node: node[1], reverse=True)

    #Create tuples of form (node, degree, color)
    color_tuples = []
    i=0
    for degree_tuple in degree_tuples:
        color_tuples.append([degree_tuple[0],degree_tuple[1],-1])
        i=i+1

    #use heuristic attempt algorithm to solve
    currentColor = 0
    color_tuples[0+offset][2]=0
    while len([x for [x,y,z] in color_tuples if z ==-1])>0:

        for node in color_tuples:
            node_index = color_tuples.index([[x,y,z] for [x,y,z] in color_tuples if x ==node[0]][0])
            if(node[2]==-1):
                nodeEdges = [[x,y] for [x,y] in edges if x==node[0] or y==node[0]]
                connections = {a for [a,b] in nodeEdges if b==node[0]}|{d for [c,d] in nodeEdges if c==node[0]}
                currentColorNodes = {x for [x,y,z] in color_tuples if z==currentColor}
                if len(connections&currentColorNodes)== 0:
                    color_tuples[node_index][2]=currentColor
        currentColor += 1
        if len([x for [x,y,z] in color_tuples if z ==-1])>0:
            initial_node_index = color_tuples.index([[x,y,z] for [x,y,z] in color_tuples if z ==-1][0])
            color_tuples[initial_node_index][2]=currentColor

    return color_tuples

def brute(nodeCount, edges):

    all_tuples = []
    for i in range(0,nodeCount):
        all_tuples.append(bruteTrivial(nodeCount,edges,i))
    min_colors = min(max([z for (x,y,z) in tup]) for tup in all_tuples )

    color_tuples = [a for a in all_tuples if max([z for (x,y,z) in a])==min_colors][0]
    # prepare the solution in the specified output format
    numColors = max([z for (x,y,z) in color_tuples])+1
    solution = [z for (x,y,z) in sorted(color_tuples, key=lambda node: node[0])]
    outputData = str(numColors) + ' ' + str(0) + '\n'
    outputData += ' '.join(map(str, solution))

    return outputData

def getUpperBound(nodeCount, edges):
   #Create tuples of form (node, degree)
    degree_tuples = []
    for node in range(0,nodeCount):
        nodeEdges = [[x,y] for [x,y] in edges if x==node or y==node]
        degree = len(nodeEdges)
        degree_tuples.append([node, degree])

    #Sort by Degree
    degree_tuples=sorted(degree_tuples, key=lambda node: node[1], reverse=True)

    #Create tuples of form (node, degree, color)
    color_tuples = []
    i=0
    for degree_tuple in degree_tuples:
        color_tuples.append([degree_tuple[0],degree_tuple[1],-1])
        i=i+1

    #use heuristic attempt algorithm to solve
    currentColor = 0
    #color_tuples[0][2]=0
    while len([x for [x,y,z] in color_tuples if z ==-1])>0:
        initial_node_index = color_tuples.index([[x,y,z] for [x,y,z] in color_tuples if z ==-1][0])
        color_tuples[initial_node_index][2]=currentColor
        if len([x for [x,y,z] in color_tuples if z ==-1])>0:
            for node in color_tuples:
                node_index = color_tuples.index([[x,y,z] for [x,y,z] in color_tuples if x ==node[0]][0])
                if(node[2]==-1):
                    nodeEdges = [[x,y] for [x,y] in edges if x==node[0] or y==node[0]]
                    connections = {a for [a,b] in nodeEdges if b==node[0]}|{d for [c,d] in nodeEdges if c==node[0]}
                    currentColorNodes = {x for [x,y,z] in color_tuples if z==currentColor}
                    if len(connections&currentColorNodes)== 0:
                        color_tuples[node_index][2]=currentColor
            currentColor += 1

    # prepare the solution in the specified output format
    return  max([z for (x,y,z) in color_tuples])+1

def firstAttempt(nodeCount, edges):
    # build a trivial solution
    # every node has its own color
    solution = range(0, nodeCount)
    solution[0] = 0
    for node in range(1,nodeCount):
        nodeColor = range(0,nodeCount)
        for edgePair in range(0,nodeCount):
            #Check if edge exists
            edge1 = (node,edgePair)
            edge2 = (edgePair,node)
            if edge1 in edges:
                nodeColor = filter(lambda a: a!= solution[edgePair], nodeColor)
            if edge2 in edges:
                nodeColor = filter(lambda a: a!= solution[edgePair], nodeColor)
        solution[node] = min(nodeColor)


    # prepare the solution in the specified output format
    outputData = str(max(solution)+1) + ' ' + str(0) + '\n'
    outputData += ' '.join(map(str, solution))

    return outputData

def secondAttempt(nodeCount, edges):
    #Create tuples of form (node, degree)
    degree_tuples = []
    for node in range(0,nodeCount):
        nodeEdges = [[x,y] for [x,y] in edges if x==node or y==node]
        degree = len(nodeEdges)
        degree_tuples.append([node, degree])

    #Sort by Degree
    degree_tuples=sorted(degree_tuples, key=lambda node: node[1], reverse=True)

    #Create tuples of form (node, degree, color)
    color_tuples = []
    i=0
    for degree_tuple in degree_tuples:
        color_tuples.append([degree_tuple[0],degree_tuple[1],-1])
        i=i+1

    #use 1st attempt algorithm to solve
    for colorTuple in color_tuples:
        nodeColor = range(0,nodeCount)
        for edgePair in range(0,nodeCount):
            #Check if edge exists
            edge1 = (colorTuple[0],edgePair)
            edge2 = (edgePair,colorTuple[0])
            vertexColor = [z for (x,y,z) in color_tuples if x==edgePair ][0]
            if edge1 in edges:
                nodeColor = filter(lambda a: a!= vertexColor, nodeColor)
            if edge2 in edges:
                nodeColor = filter(lambda a: a!= vertexColor, nodeColor)
        colorTuple[2] = min(nodeColor)

    # prepare the solution in the specified output format
    numColors = max([z for (x,y,z) in color_tuples])+1
    solution = [z for (x,y,z) in sorted(color_tuples, key=lambda node: node[0])]
    outputData = str(numColors) + ' ' + str(0) + '\n'
    outputData += ' '.join(map(str, solution))

    return outputData

def thirdAttempt(nodeCount, edges, offset):
    #Create tuples of form (node, degree)
    color_tuples = []
    neighbors = []
    for node in range(0,nodeCount):
        emptySet = set()
        neighbors.append(emptySet)
    for edge in edges:
        neighbors[edge[0]].add(edge[1])
        neighbors[edge[1]].add(edge[0])
    for node in range(0,nodeCount):
        color_tuples.append([node, len(neighbors[node]), -1])

    #Sort by Degree
    color_tuples=sorted(color_tuples, key=lambda node: node[1], reverse=True)

    #use heuristic attempt algorithm to solve
    currentColor = 0
    #color_tuples[0][2]=0
    hasRun = False
    while len([x for [x,y,z] in color_tuples if z ==-1])>0:
        startNode = 0
        startNode = 0 if hasRun else offset
        initial_node_index = color_tuples.index([[x,y,z] for [x,y,z] in color_tuples if z ==-1][startNode])
        color_tuples[initial_node_index][2]=currentColor
        initialNode = color_tuples[initial_node_index][0]
        hasRun=True
        if len([x for [x,y,z] in color_tuples if z ==-1])>0:
            #find list of non-neighbors to initial node
            NonNeighbors = range(0,nodeCount)
            for node in [[x,y,z] for [x,y,z] in color_tuples if z !=-1]:
                NonNeighbors.remove(node[0])
            ColorNeighbors=set()
            #NonNeighbors.remove(initialNode)
            for i in neighbors[initialNode]:
                if i in NonNeighbors: NonNeighbors.remove(i)
                ColorNeighbors.add(i)
            while len(NonNeighbors)>0:
                #find the non-neighbor with the most neighbors in common with initial node
                biggestNN=-1
                biggestSize=-1
                for NN in NonNeighbors:
                    if len(ColorNeighbors & neighbors[NN]) > biggestSize:
                        biggestNN=NN
                        biggestSize = len(ColorNeighbors & neighbors[NN])

                #color this node same as initial node
                node_index = color_tuples.index([[x,y,z] for [x,y,z] in color_tuples if x ==biggestNN][0])
                color_tuples[node_index][2] = currentColor
                NonNeighbors.remove(biggestNN)
                for i in neighbors[biggestNN]:
                    if i in NonNeighbors: NonNeighbors.remove(i)
                    ColorNeighbors.add(i)
            currentColor += 1

    # prepare the solution in the specified output format
    numColors = max([z for (x,y,z) in color_tuples])+1
    solution = [z for (x,y,z) in sorted(color_tuples, key=lambda node: node[0])]
    outputData = str(numColors) + ' ' + str(0) + '\n'
    outputData += ' '.join(map(str, solution))

    return outputData

def thirdRecursive(nodeCount, edges):
    bestSize = nodeCount
    bestOut = ''

    for i in range(0,10):
        tempOut = thirdAttempt(nodeCount, edges, i)
        outLines=tempOut.splitlines()
        size = int(outLines[0].split()[0])
        if size<bestSize:
            bestSize=size
            bestOut=tempOut
    return bestOut


def solverAttempt (nodeCount, edges):
    #Create solver model file
    upperBound=getUpperBound(nodeCount,edges)
    path='C:\\Solver\\bin\\'
    filename='test_'+str(nodeCount)+'_'+str(datetime.datetime.now().date())+'_'+ str(datetime.datetime.now().hour)+'_'+str(datetime.datetime.now().minute)+'_'+str(datetime.datetime.now().second)+'.mzn'
    solverFile = open(path+filename,'w')
    solverFile.write('%G12 Solver File for '+str(nodeCount)+' nodes run at '+ str(datetime.datetime.now())+'\n')
    solverFile.write('int: nc='+str(nodeCount)+'; \n')
    solverFile.write('\n')
    solverFile.write('array[0..nc] of var 0..'+str(upperBound)+': color; \n')
    solverFile.write('\n')

    solverFile.write('constraint color[0]=0;')

    solverFile.write('\n')

    for edge in edges:
        solverFile.write('constraint color['+str(edge[0])+']!=color['+str(edge[1])+']; \n')

    solverFile.write('\n')
    solverFile.write('solve minimize max(color); \n')
    solverFile.write('')
    solverFile.write('output [ show(color[i]) ++"\\n" | i in 0..nc]; \n')
    solverFile.close()

    #run solver and get solution
    import os
    from subprocess import check_output as qx

    os.chdir('c:\\Solver\\bin')
    cmd = r'mzn-g12fd.bat '+filename
    output = qx(cmd)
    os.chdir('C:\\Users\\keseaman\\Dropbox\\Coursera\\Discrete Optimization\\coloring')
    lines=output.splitlines();

    #format output and return
    numColors = int(max(lines [:nodeCount]))+1
    opt = 0
    if lines [-1]=='==========': opt = 1

    outputData = str(numColors) + ' ' + str(opt) + '\n'
    outputData += ' '.join(map(str, lines[:nodeCount]))
    return outputData;

def solverAttempt2 (nodeCount, edges):


    try:
        #Create a new Model
        m = Model("myModel")
    except GurobiError:
            print 'Error Reported'

    return 'Test'

def solveIt(inputData):
    # Modify this code to run your optimization algorithm

    # parse the input
    lines = inputData.split('\n')

    firstLine = lines[0].split()
    nodeCount = int(firstLine[0])
    edgeCount = int(firstLine[1])
    #
    #
    #
    edges = []
    for i in range(1, edgeCount + 1):
        line = lines[i]
        parts = line.split()
        edges.append((int(parts[0]), int(parts[1])))

    #outputData = trivialSolution(nodeCount)
    #outputData = firstAttempt(nodeCount, edges)
    #outputData = secondAttempt(nodeCount, edges)
    #outputData = secondTrivial(nodeCount, edges)
    #outputData = brute(nodeCount, edges)
    #outputData = solverAttempt(nodeCount, edges)
    #outputData = solverAttempt2(nodeCount, edges)
    #outputData = thirdAttempt(nodeCount, edges, 0)
    outputData = thirdRecursive(nodeCount, edges)

    return outputData



import sys

if __name__ == '__main__':
    if len(sys.argv) > 1:
        fileLocation = sys.argv[1].strip()
        inputDataFile = open(fileLocation, 'r')
        inputData = ''.join(inputDataFile.readlines())
        inputDataFile.close()
        print solveIt(inputData)
    else:
        print 'This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/gc_4_1)'

