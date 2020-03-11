#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np

readsFile = open("reads.fq", "r")
readsLines = readsFile.readlines()
readsFile.close()

reads = []
for i in range(1, len(readsLines), 4):
    
    reads.append(readsLines[i][:-1])

#read1 = "ATTCGTCAGTATGGATCAT"
#read2 = "ATGATCATCCATGTC"

#read1 = "GTCCATAGATCC"
#read2 = "CAAGACCACTTGAG"

#read1 = reads[0]
#read2 = reads[1]


def chevauchement(seq1, seq2):
    read1 = seq1
    read2 = seq2
    
    matrix = np.ones((len(read2)+1, len(read1)+1))
    matrix[:,0] = 0
    matrix[0] = 0

    chemins = np.zeros((len(read2)+1, len(read1)+1))
    chemins = chemins - 1

    for i in range(1, len(read2)+1):
        for j in range(1, len(read1)+1):
            diagoValue = -4
            if(read1[j-1] == read2[i-1]):
                diagoValue = 4

            matrix[i,j] = max(matrix[i-1,j]-8, matrix[i, j-1]-8, matrix[i-1, j-1] + diagoValue)

            if(matrix[i,j] == matrix[i-1,j]+8):
                chemins[i,j] = 0
            elif(matrix[i,j] == matrix[i,j-1]+8):
                chemins[i,j] = 1
            else:
                chemins[i,j] = 2
    #print(matrix)
    #print("\n")

    #print(chemins)
    #print("\n")

    lastColumn = list(matrix[:,len(read1)])

    #print(lastColumn)
        
    response = traceback(np.argmax(lastColumn), matrix, chemins, read1, read2)
    return (max(lastColumn), response, len(response[0]))



def traceback(startIndex, matrix, chemins, read1, read2):
    tracebackList1 = []
    tracebackList2 = []
    
    i = startIndex
    j = len(read1)
    while(chemins[i,j]!=-1):
        if(chemins[i,j]==0):
            tracebackList1.append("-")
            tracebackList2.append(read2[i-1])
            i = i-1
        elif(chemins[i,j]==1):
            #tracebackList.append("D")
            tracebackList1.append(read1[j-1])
            tracebackList2.append("-")
            j = j-1
        else:
            #tracebackList.append("M")
            #print(matrix)
            #print(i,j)
            #print(read2)
            tracebackList1.append(read1[j-1])
            tracebackList2.append(read2[i-1])
            i = i-1
            j = j-1
            
    tracebackList1.reverse()
    tracebackList2.reverse()
    
    #traceback(localMin(lastColumn, len(lastColumn))-1)
    #traceback(np.argmax(lastColumn))
    #print(tracebackList1)
    #print(tracebackList2)
    
    return(tracebackList1, tracebackList2)

def alignAll():
    alignments = np.zeros((len(reads), len(reads)))
    alignments = alignments -1
    
    counter = 0
    
    for i in range(len(reads)):
        for j in range(len(reads)):
            if(i != j):
                response = chevauchement(reads[i], reads[j])
                alignments[i,j] = int(response[0])
            counter += 1
            print(counter)
                
    return alignments

alignmentMatrix = alignAll()
print(alignmentMatrix)



# In[22]:


print(len(reads[0]))

outputFile = open("output.txt", "w")

line = ""

for i in range(len(alignmentMatrix[0])):
    for j in range(len(alignmentMatrix[0])):
        if(i != j):
            line += "seq" + str(i+1) + "\tseq" + str(j+1) + "\t" + str(alignmentMatrix[i,j]) + "\n"
outputFile.write(line)


# In[21]:


outputFileTresh = open("outputTresh.txt", "w")

line = ""

for i in range(len(alignmentMatrix[0])):
    for j in range(len(alignmentMatrix[0])):
        if(i != j and int(alignmentMatrix[i,j]) >= 80):
            line += "seq" + str(i+1) + "\tseq" + str(j+1) + "\t" + str(alignmentMatrix[i,j]) + "\n"
outputFileTresh.write(line)

