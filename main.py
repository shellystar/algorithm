
ReadsList = []

ReadsLen = 100
#ReadsLen = 8

Graph = {}

KLEN = 19
#KLEN = 5

STARTKEY = []
CHECKKEY = []

TIPSLENLIMIT = ReadsLen - KLEN

CONTIGS = []

DATA4CNT = 0

class diyError(Exception):
    def __init__(self, errorinfo):
        super().__init__(self)
        self.errorinfo = errorinfo

class KMERNODE():
    def __init__(self):
        self.outDegree = 0
        self.inDegree = 0
        self.coverage = 1
        self.next = {}  # elements are dictionary consisting of next node's key as key and weights of this edge as value
        self.last = {}  # same structure as next
        ## store which contig it belongs to, identical to the index of start nodes, which is able to reach this node, in STAETKEY
        #self.configNum = []
    
    def addOutEdge(self, nextKey, edgeWeight):
        if nextKey in self.next.keys():
            self.next[nextKey] += edgeWeight
            return
        self.next[nextKey] = edgeWeight
        self.outDegree += 1
    
    def addInEdge(self, lastKey, edgeWeight):
        if lastKey in self.last.keys():
            self.last[lastKey] += edgeWeight
            return
        self.last[lastKey] = edgeWeight
        self.inDegree += 1

    def delOutEdge(self, nextKey):
        if nextKey in self.next.keys():
            del(self.next[nextKey])
            self.outDegree -= 1
            return 0
        else:
            return -1
    
    def delInEdge(self, lastKey):
        if lastKey in self.last.keys():
            del(self.last[lastKey])
            self.inDegree -= 1
            return 0
        else:
            return -1

    def updateEdge(self, oldKey, newKey, flag):
        """
        flag 0 refers to update oldKey in last nodes
        flag 1 refers to update oldKey in next nodes
        """
        if flag == 0 and oldKey in self.last.keys():
            weight = self.last.pop(oldKey)
            self.last[newKey] = weight
            return 0
        elif flag == 1 and oldKey in self.next.keys():
            weight = self.next.pop(oldKey)
            self.next[newKey] = weight
            return 0
        else:
            return -1

def complement(originStr):
    tmp = list(originStr)
    for i in range(len(tmp)):
        if tmp[i] == 'A':
            tmp[i] = 'T'
        elif tmp[i] == 'T':
            tmp[i] = 'A'
        elif tmp[i] == 'C':
            tmp[i] = 'G'
        elif tmp[i] == 'G':
            tmp[i] = 'C'
    return ''.join(tmp)

def reverse(originStr):
    originStr = originStr[-1::-1]
    return originStr

def readsOrientation(reads):
    global Graph

    hitTimes = [0] * len(reads)
    for i in range(len(reads)):
        for k in range(ReadsLen - KLEN + 1):
            if reads[i][k : k + KLEN] in Graph.keys():
                hitTimes[i] += 1
    return hitTimes.index(max(hitTimes))

def getNode(key):
    global Graph

    if key in Graph.keys():
        Graph[key].coverage += 1
        return Graph[key]
    else:
        Graph[key] = KMERNODE()
        return Graph[key]

def loadData(fastaFile):
    with open(fastaFile, 'r') as f:
        cnt = 1
        for line in f:
            if cnt % 2 == 0:
                ReadsList.append(line.split('\n')[0])
            cnt += 1

def saveData(fastaFile):
    global CONTIGS

    f = open(fastaFile, 'w')
    cnt = 1
    N = len(CONTIGS)
    print ("Total contigs: ", N)
    for i in range(N):
        if len(CONTIGS[i]) < 100:
            continue
        f.write(">contig_" + str(cnt) + "/" + str(N) + "\n")
        f.write(str(CONTIGS[i]) + '\n')
        cnt += 1
    f.close()

def constructGraph():
    global Graph, KLEN, ReadsList, ReadsLen

    N = len(ReadsList)
    for i in range(len(ReadsList)):
        reads = ['reads'] * 4
        reads[0] = ReadsList[i]
        reads[1] = complement(reads[0])
        reads[2] = reverse(reads[0])
        reads[3] = complement(reads[2])

        # choose one orientation which has the biggest hit times and construct nodes and edges
        for j in range(4):
            read = reads[j]
            for k in range(ReadsLen - KLEN + 1):
                kmer = read[k : k + KLEN]
                node = getNode(kmer)
                if k != 0:
                    lastKmer = read[k-1 : k-1 + KLEN]
                    lastNode = Graph[lastKmer]
                    lastNode.addOutEdge(kmer, 1)
                    node.addInEdge(lastKmer, 1)
            print ("Constructing graph. ", i, "/", N, " complete")

    del(ReadsList)

def reduceLowCoverage(threshold):
    global Graph

    for key in Graph:
        if Graph[key].coverage <= threshold:
            for k in Graph[key].last:
                Graph[k].delOutEdge(key)
            for k in Graph[key].next:
                Graph[k].delInEdge(key)
            Graph[key] = 0

def saveGraph(graphFile):
    global Graph
    with open(graphFile, 'w') as f:
        for key in Graph:
            if Graph[key] == 0:
                continue
            node = Graph[key]
            f.write(str(key) + ' ' + str(node.coverage) + '\n')
            for k in node.last:
                f.write(str(k) + ',' + str(node.last[k]) + ' ')
            f.write('\n')
            for k in node.next:
                f.write(str(k) + ',' + str(node.next[k]) + ' ')
            f.write('\n')

def loadGraph(graphFile):
    global Graph

    with open(graphFile, 'r') as f:
        cnt = 1
        for line in f:
            if cnt % 3 == 1:
                line = line.split('\n')[0].split(' ')
                nodeKey = line[0]
                Graph[nodeKey] = KMERNODE()
                Graph[nodeKey].coverage = int(line[1])
            elif cnt % 3 == 2:
                line = line.split(' ')
                for i in range(len(line)-1):
                    Graph[nodeKey].addInEdge(line[i].split(',')[0], int(line[i].split(',')[1]))
            elif cnt % 3 == 0:
                line = line.split(' ')
                for i in range(len(line)-1):
                    Graph[nodeKey].addOutEdge(line[i].split(',')[0], int(line[i].split(',')[1]))
            cnt += 1

def currentTerminal(flag):
    """
    flag indicates which kind of terminals to find
    1 for start point(inDegree is 0)
    0 for end point(outDegree is 0)
    """
    keyList = []

    if flag == 1:
        for key in Graph:
            if Graph[key] == 0:
                continue
            if Graph[key].inDegree == 0:
                keyList.append(key)
    else:
        for key in Graph:
            if Graph[key] == 0:
                continue
            if Graph[key].outDegree == 0:
                keyList.append(key)
    return keyList

def mergeNext(thisKey, nextKey):
    global CHECKKEY

    thisNode = Graph[thisKey]
    nextNode = Graph[nextKey]

    # judge whether be able to merge or not
    if nextKey not in thisNode.next: ##or thisNode.outDegree > 1:
        return -1
    
    # merge KMERNODE
    thisNewKey = thisKey + nextKey[KLEN-1 : ]
    thisNewNode = Graph.pop(thisKey)
    Graph[thisNewKey] = thisNewNode
    # update nodes which points to thisNewNode
    for key in thisNewNode.last:
        lastNode = Graph[key]
        lastNode.updateEdge(thisKey, thisNewKey, 1)

    # update thisnewnode's out edges and coverage information, and nextnode's in edges information
    thisNewNode.delOutEdge(nextKey)
    nextNode.delInEdge(thisKey)
    for key in nextNode.next:
        thisNewNode.addOutEdge(key, nextNode.next[key])
        Graph[key].updateEdge(nextKey, thisNewKey, 0)
    #thisNode.coverage = min(thisNode.coverage, nextNode.coverage)

    # check whether this key is in the CHECKKEY list, if it is, change the key value in CHECKKEY list
    if thisKey in CHECKKEY:
        idx = CHECKKEY.index(thisKey)
        CHECKKEY[idx] = thisNewKey
    if thisKey in STARTKEY:
        idx = STARTKEY.index(thisKey)
        STARTKEY[idx] = thisNewKey
    
    # if nextnode's indegree is 0, delete it in the Graph
    #if nextNode.inDegree == 0:
    #    del(Graph[nextKey])
    return thisNewKey

def solveCircle(startKey, nextKey):
    """
    dissolve circles
    merge all nodes appearing in circles as one node to dissolve it
    """
    global Graph

    # step1: merge all nodes to node B except startNode A
    # step2: merge B with A --> BA, reupdate all in edges of out edges of BA, and delete all out edges of BA
    # step3: merge A with BA

    # step1
    bPtrKey = nextKey
    fPtrKey = list(Graph[bPtrKey].next.keys())[0]
    while(fPtrKey != startKey):
        bPtrKey = mergeNext(bPtrKey, fPtrKey)
        fPtrKey = list(Graph[bPtrKey].next.keys())[0]

    # step2
    bPtrKey = mergeNext(bPtrKey, fPtrKey)
    bPtrNode = Graph[bPtrKey]
    nextKeysList = list(bPtrNode.next.keys())
    for k in nextKeysList:
        Graph[k].updateEdge(bPtrKey, startKey, 0)
        bPtrNode.delOutEdge(k)

    # step3
    circleKey = mergeNext(startKey, bPtrKey)

    if Graph[circleKey].outDegree == 1 and Graph[circleKey].inDegree == 1:
        return 0
    else:
        raise diyError('SolveCircle Error')
    
    print ("call solveCircle.")
    
def solveBubble(startKey, endKey):
    """
    dissolve bubbles
    delete the branch with lower edge weight
    """
    global Graph

    startNode = Graph[startKey]
    lightestEdgeKey = min(startNode.next.items(), key=lambda x: x[1])[0]
    
    # tour along lightestEdgeKey and delete this branch
    delKey = lightestEdgeKey
    startNode.delOutEdge(delKey)
    Graph[delKey].delInEdge(startKey)
    
    bKey = startKey
    while(delKey != endKey):
        delNode = Graph[delKey]
        bKey = delKey
        delKey = list(delNode.next.keys())[0]
    
    Graph[bKey].delOutEdge(delKey)
    Graph[delKey].delInEdge(bKey)
    print ("call solveBubble.")

def solveTips(juncKey, tipsKey):
    """
    delete tips branch
    """
    global Graph

    # tour along tipsKey and delete tips branch
    """delKey = juncKey
    nextKey = tipsKey
    while(Graph[delKey].outDegree != 0):
        delNode = Graph[delKey]
        nextNode = Graph[nextKey]
        delNode.delOutEdge(nextKey)
        nextNode.delInEdge(delKey)
        
        #if delNode.inDegree == 0 and delNode.outDegree == 0:
        #    del (Graph[delKey])
        delKey = nextKey
        if nextNode.outDegree:
            nextKey = list(nextNode.next.keys())[0]
    
    delNode = Graph[delKey]
    #if delNode.inDegree == 0 and delNode.outDegree == 0:
    #    del Graph[delKey]"""

    juncNode = Graph[juncKey]
    tipsNode = Graph[tipsKey]
    juncNode.delOutEdge(tipsKey)
    tipsNode.delInEdge(juncKey)

    print ("call solveTips.")

def solveCross(crossKey):
    """
    match branches with weights of edges
    merge to dispose of cross nodes
    """
    global Graph, CROSSCALL

    #if CROSSCALL == 584:
    #    print ("here comes!")

    # merge cross nodes with every nodes pointing to it, delete unmatched next nodes(match by weights of edges)
    """crossNode = Graph[crossKey]
    lastKeysList = list(crossNode.last.keys())
    for key in lastKeysList:
        branchWeight = crossNode.last[key]

        print ("cross solve: ", crossKey, " ", key)

        key = mergeNext(key, crossKey)
        branchNode = Graph[key]
        branchNextKeys = list(branchNode.next.keys())
        for out in branchNextKeys:
            if branchNode.next[out] != branchWeight:
                branchNode.delOutEdge(out)
                Graph[out].updateEdge(key, crossKey, 0)  # not delete out edges' in edge
            else:
                crossNode.delOutEdge(out)"""
    crossNode = Graph[crossKey]
    crossLastKeys = list(crossNode.last.keys())
    crossNextKeys = list(crossNode.next.keys())
    lastWeights = []
    nextWeights = []

    for key in crossLastKeys:
        lastWeights.append(crossNode.last[key])

    for key in crossNextKeys:
        nextWeights.append(crossNode.next[key])

    for key in crossLastKeys:
        newKey = mergeNext(key, crossKey)
        crossLastKeys[crossLastKeys.index(key)] = newKey
        for k in crossNextKeys:
            Graph[k].updateEdge(newKey, crossKey, 0)
    
    branchNum = min(len(crossLastKeys), len(crossNextKeys))
    for i in range(branchNum):
        key = crossLastKeys[i]
        
        weightsDist = []
        for j in range(branchNum):
            weightsDist.append(abs(lastWeights[i] - nextWeights[j]))
        minDistIdx = weightsDist.index(min(weightsDist))
        nextWeights[minDistIdx] = -10000

        matchKey = crossNextKeys[minDistIdx]
        Graph[matchKey].updateEdge(crossKey, key, 0)

        for k in crossNextKeys:
            if k != matchKey:
                Graph[key].delOutEdge(k)    

    print ("call solveCross")      
    #CROSSCALL += 1    

def junctionNodes_back(key):
    global Graph, CHECKKEY

    node = Graph[key]
    # in degree > 1: cross, circle 
    # out degree > 1: tips, multi contigs and bubbles
    if node.inDegree > 1:
        # By merging all nodes in route, put all probable branches under out degree > 1 situation
        while(node.outDegree == 1):
            key = mergeNext(key, list(node.next.keys())[0])
            if key == -1:
                raise diyError("Merge error!")
            node = Graph[key]
    
    if node.outDegree > 1:
        # tour each branch
        routeCnt = []   # store the length of each branch route
        branchKey = []
        bubbleEndKey = 0
        nodeNextList = list(node.next.keys())
        for juncNextKey in nodeNextList:
            branchKey.append(juncNextKey)
            juncNextNode = Graph[juncNextKey]
            cnt = 1
            tourKey = juncNextKey
            tourNode = juncNextNode
            while(tourNode.outDegree == 1 and tourNode.inDegree == 1):
                tourKey = list(tourNode.next.keys())[0]
                tourNode = Graph[tourKey]
                cnt += 1
            routeCnt.append(cnt)

            # terminate condition: indegree != 1, only two circumstances: circle and bubble
            if tourNode.inDegree > 1:
                # circle
                if tourNode == node:
                    solveCircle(key, juncNextKey)
                    routeCnt[-1] = -1   # mark routeCnt as -1, indicating that this branch has been solved
                    return
                # bubble
                else:
                    if not bubbleEndKey:
                        bubbleEndKey = tourKey
                    elif bubbleEndKey == tourKey:
                        solveBubble(key, tourKey)
                        routeCnt[-1] = -1
                        return
            
        shortestLen = min(routeCnt)
        if shortestLen > 0:
            # terminate condition: outdegree == 0, possible circumstances: cross point, tips, multi contigs junction nodes
            # cross point
            if node.inDegree > 1 and node.inDegree == node.outDegree:
                solveCross(key)
                return
            # tips
            elif shortestLen < TIPSLENLIMIT:
                tipsKey = branchKey[routeCnt.index(shortestLen)]
                solveTips(key, tipsKey)  
                return
            # multi contigs
            else:
                for k in node.next:
                    CHECKKEY.append(k)
                return

def junctionNodes(key):
    """
    only solve cross nodes and circle
    """
    global Graph, CHECKKEY

    node = Graph[key]
    # in degree > 1: cross, circle or bubble end
    # only solve circle and cross 
    if node.inDegree > 1:
        # By merging all nodes in route, put all probable branches under out degree > 1 situation
        while(node.outDegree == 1):
            key = mergeNext(key, list(node.next.keys())[0])
            if key == -1:
                raise diyError("Merge error!")
            node = Graph[key]
        
        if node.outDegree <= 1:
            return
        
        # tour each branch
        nodeNextList = list(node.next.keys())
        for juncNextKey in nodeNextList:
            juncNextNode = Graph[juncNextKey]
            tourKey = juncNextKey
            tourNode = juncNextNode
            while(tourNode.outDegree == 1 and tourNode.inDegree == 1):
                tourKey = list(tourNode.next.keys())[0]
                tourNode = Graph[tourKey]

            # terminate condition: indegree != 1, circle
            if tourNode.inDegree > 1:
                # circle
                if tourNode == node:
                    solveCircle(key, juncNextKey)
                    return
        
        # already tour each branch and no circle nodes, cross
        solveCross(key)
        return
    if node.outDegree > 1:
        return

def fixGraph():
    global Graph, STARTKEY, CHECKKEY

    STARTKEY = currentTerminal(1)
    
    CHECKKEY = STARTKEY[:]

    return
    checkPtNum = len(CHECKKEY)
    print ("origin checknum: ", checkPtNum)
    presentNum = 0

    
    while(presentNum != checkPtNum):
        for i in range(presentNum, checkPtNum):
            bPtrKey = -1
            fPtrKey = CHECKKEY[i]
            fPtrNode = Graph[fPtrKey]
            while(fPtrNode.outDegree != 0):
                cnt = 0
                if fPtrNode.outDegree > 1 or fPtrNode.inDegree > 1:
                    #junctionNodes_back(fPtrKey)
                    junctionNodes(fPtrKey)
                    if fPtrKey in Graph:
                        if Graph[fPtrKey].outDegree > 1:
                            break
                    bPtrKey = -1
                    try:
                        fPtrKey = CHECKKEY[i]
                        fPtrNode = Graph[fPtrKey]
                    except KeyError:
                        print ("-----------------------start point {:}, ------------------------")
                else:
                    bPtrKey = fPtrKey
                    fPtrKey = list(fPtrNode.next.keys())[0]
                    #if fPtrKey == 'CCATTCGCATAGCGGGAGC':
                    #    print ("error here")
                    fPtrNode = Graph[fPtrKey]
            presentNum += 1
        checkPtNum = len(CHECKKEY)
        print ("update checkptnum: ", checkPtNum)

def tour(key, contig):
    global Graph, KLEN, CONTIGS, DATA4CNT

    #if DATA4CNT >= 10:
    #    return

    stack = []
    visited = {}
    stack.append((key, contig))
    # tour until outdegree > 1: push into stack
    # outdegree == 0: add an contig, pop a stack
    while(len(stack) != 0):
        key, contig = stack.pop()
        if key in visited:
            continue
        visited[key] = True
        node = Graph[key]
        if node.outDegree == 1:
            contig += key[KLEN-1 : ]
            key = list(node.next.keys())[0]
            stack.append((key, contig))
        elif node.outDegree == 0:
            contig += key[KLEN-1 : ]
            #CONTIGS.append(contig)
            print (len(contig))
            if len(contig) < 100000 and len(contig) > 5000:
                fw = open('./contig/data4/version3.0Contig_1_10klen.fasta', 'a')
                fw.write(">contig\n")
                fw.write(contig + '\n')        
                fw.close()        
        else:
            contig += key[KLEN-1 : ]
            for k in node.next:
                stack.append((k, contig))
    return

def getContigs():
    global STARTKEY, Graph, CONTIGS
    for i in range(len(STARTKEY)):
        key = STARTKEY[i]
        contig = key[: KLEN-1]
        tour(key, contig)

def data4GetContigs():
    global STARTKEY, Graph, CONTIGS
    for i in range(len(STARTKEY)):
        key = STARTKEY[i]
        contig = key[:KLEN-1]
        DFSonce(key, contig)

def DFSonce(key, contig):
    global Graph, CONTIGS

    while(Graph[key].outDegree != 0 and len(contig) < 99999):
        contig += key[KLEN-1: ]
        if Graph[key].outDegree == 1:
            key = list(Graph[key].next.keys())[0]
        else:
            keyList = list(Graph[key].next.keys())
            weightList = [0] * Graph[key].outDegree
            for i in range(len(keyList)):
                weightList[i] = Graph[key].next[keyList[i]]
            maxWidx = weightList.index(max(weightList))
            key = keyList[maxWidx]
    contig += key[KLEN-1: ] 
    if len(contig) > 50000 and len(contig) < 100000:
        CONTIGS.append(contig)
    if len(contig) >= 100000:
        CONTIGS.append(contig[:99999])
        if len(contig) - 99999 > 50000:
            CONTIGS.append(contig[9999:])
    print (len(contig))

if __name__ == '__main__':
    #loadData('./data/data4/short_1.fasta')
    #loadData('./data/data4/short_2.fasta')
    #loadData('./testMultiContig_reads.fasta')

    #constructGraph()
    #reduceLowCoverage(1)
    #saveGraph('./graph/data4/graphData_version3.0_19klen_data2.txt')
    loadGraph('./graph/data4/graphData_version3.0_19klen_data4.txt')
    print (len(Graph))
    fixGraph()
    #getContigs()
    data4GetContigs()

    #saveData('./testMultiContig_contig.fasta')
    saveData('./contig/data4/version3.0Contig_1_19klen.fasta')
