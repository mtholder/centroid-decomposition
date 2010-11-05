#include "ncl/nxsmultiformat.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <list>


void printHelp(std::ostream & out);
bool treeReadCallback(NxsFullTreeDescription &, void *, NxsTreesBlock *);

// some globals
std::string gVersionString("centroid_decomp version 0.0.1");
std::string gErrMsgPrefix("centroid_decomp: ");
MultiFormatReader * gNexusReader = 0L;
bool gTaxaBlockWritten = false;
int gCurrTreeIndex = 0;
NxsString gErrorMessage;
NxsSimpleNode bogusNode(0L, 0.0);
unsigned gTaxaBlockNumLeaves = 0;
NxsTaxaBlock * gPrevTaxaBlock = 0L;
NxsSimpleNode * gRoot = 0L;
long gNumNodes = 0;
long gNumLeaves = 0;
bool gDebugging = false;
enum TreeSweepDirection {
    LEFT_DIR_BIT = 1,
    RIGHT_DIR_BIT = 2,
    LEFT_OR_RIGHT = 3,
    BELOW_DIR_BIT = 4,
    LEFT_OR_BELOW = 5,
    RIGHT_OR_BELOW = 6,
    ALL_DIR_BITS = 7,
    THIS_NODE_DIR_BIT = 8
};
enum TreeDirection {
    LEFT_DIR = 1,
    RIGHT_DIR = 2,
    BELOW_DIR = 4,
    THIS_NODE_DIR = 8
};

class LeafPathElement {
    public:
    LeafPathElement(const NxsSimpleNode *nd, TreeDirection d) {
        assert(nd);
        assert(nd->IsTip());
        leaf = nd;
        iScore = 1;
        dirToNext = d;
        indexInNext = -1;
    }

    LeafPathElement(const LeafPathElement &nextEl, TreeDirection d, unsigned index) {
        leaf = nextEl.leaf;
        dirToNext = d;
        assert(index < INT_MAX);
        indexInNext = (int)index;
    }
    
    bool operator<(const LeafPathElement & other) const {
        return this->iScore < other.iScore;
    }
    public: // should be private
        const NxsSimpleNode * leaf;
        long iScore;
        TreeDirection dirToNext;
        int indexInNext;
};
typedef std::vector<LeafPathElement> LPECollection;
typedef LPECollection::const_iterator LPECollectionConstIt;
typedef LPECollection::iterator LPECollectionIt;
std::ostream & operator<<(std::ostream &, const LeafPathElement &);

class EdgeDecompInfo {
    public:
        LPECollection closestLeavesAbove;
        LPECollection closestLeavesBelow;
};

class NdBlob {
    public:
    
        NdBlob(bool parentsLeftC): isParentsLeftChild(parentsLeftC) {
            Reset();
        }
        
        long GetNumLeavesBelow() const {
            return gNumNodes - this->numLeavesAboveEdge;
        }
        void Reset();
        void pop();
        
        long numLeavesAboveEdge;
        int ndDirWRTParent;

        EdgeDecompInfo fullEdgeInfo;
        long numActiveLeavesAboveEdge;
        TreeSweepDirection activeLeafDir;
        std::stack<EdgeDecompInfo*> edgeInfoStack;
        std::stack<TreeSweepDirection> activeLeafDirStack;
        EdgeDecompInfo * activeEdgeInfo;
        
        LPECollection lpeScratch1;
        LPECollection lpeScratch2;
        
        bool isParentsLeftChild;
        
//        long numLeavesBelowEdge;
};



inline void NdBlob::Reset() {
    this->numLeavesAboveEdge = -1;
    assert(edgeInfoStack.size() < 2);
    if (edgeInfoStack.size() == 1) {
        delete edgeInfoStack.top();
        edgeInfoStack.pop();
    }
    assert(edgeInfoStack.empty());
    activeEdgeInfo = &(this->fullEdgeInfo);
    //this->numLeavesBelowEdge = -1;
}

inline void NdBlob::pop() {
    if (activeEdgeInfo != &(this->fullEdgeInfo)) {
        delete activeEdgeInfo;
    }
    edgeInfoStack.pop();
    activeEdgeInfo = edgeInfoStack.top();
    activeLeafDirStack.pop();
    activeLeafDir  = activeLeafDirStack.top();
}



inline std::ostream & operator<<(std::ostream & o, const TreeDirection &d) {
    if (d == LEFT_DIR)
        o << "LEFT_DIR";
    else if (d == RIGHT_DIR)
        o << "RIGHT_DIR";
    else if (d == BELOW_DIR)
        o << "BELOW_DIR";
    else 
        o << "THIS_NODE_DIR";
    return o;
}

inline std::ostream & operator<<(std::ostream & o, const LeafPathElement &lpe) {
    o << "LeafPathElement(leaf= "<< (long)lpe.leaf << ", dirToNext= " << lpe.dirToNext << ", indexInNext= " << lpe.indexInNext << ")";
    return o;
}


unsigned gMaxNumLeafPaths = 0;

// user-controlled options as global variables...
bool gVerbose = false;
long gLeafSetIntersectionSize = 50; // we demand that the leaf set intersection size be exactly THIS_NODE_DIR
long gFloorHalfIntersectionSize = 25; // we shoot for this many above and below the centroid

long gMaxSubProblemSize = 250;
bool gUseEdgeLengths = false;
bool gIntercalate = false;

std::ostream * gOutputStream = &(std::cout);

inline bool IsRootChild(const NxsSimpleNode &nd) {
    const NxsSimpleNode *lc = nd.GetParent();
    if (lc == 0L)
        return false;
    return (lc->GetParent() == 0L);
}

inline const NxsSimpleNode * LeftChild(const NxsSimpleNode &nd) {
    const NxsSimpleNode *lc = nd.GetFirstChild();
    assert(lc);
    return lc;
}

inline const NxsSimpleNode * RightChild(const NxsSimpleNode &nd) {
    const NxsSimpleNode *lc = nd.GetFirstChild();
    assert(lc);
    const NxsSimpleNode *rc = lc->GetNextSib();
    assert(rc);
    return rc;
}

inline const NxsSimpleNode * Parent(const NxsSimpleNode &nd) {
    const NxsSimpleNode *lc = nd.GetParent();
    assert(lc);
    return lc;
}

inline NdBlob * LeftBlob(const NxsSimpleNode &nd) {
    const NxsSimpleNode *lc = LeftChild(nd);
    assert(lc);
    assert(lc->scratch);
    return (NdBlob *)lc->scratch;
}

inline NdBlob * RightBlob(const NxsSimpleNode &nd) {
    const NxsSimpleNode *lc = RightChild(nd);
    assert(lc);
    assert(lc->scratch);
    return (NdBlob *)lc->scratch;
}

inline NdBlob * ParentBlob(const NxsSimpleNode &nd) {
    const NxsSimpleNode *lc = Parent(nd);
    assert(lc);
    assert(lc->scratch);
    return (NdBlob *)lc->scratch;
}



void mergePathElementLists(LPECollection &peList,
                           TreeDirection firDir, 
                           const LPECollection & firSource,
                           TreeDirection secDir, 
                           const LPECollection & secSource
                           ) {
    //std::cerr << "mergePathElementLists insizes = " << peList.size() << ", " << firSource.size() << ", " << secSource.size()   << std::endl;
    peList.reserve(gMaxNumLeafPaths);
    assert(!gIntercalate);
    if (firSource.empty()) {
        LPECollectionConstIt secIt = secSource.begin();
        for (unsigned ind = 0; secIt != secSource.end(); ++secIt, ++ind) {
            peList.push_back(LeafPathElement(*secIt, secDir, ind));
        }
        return;
    }
    if (secSource.empty()) {
        LPECollectionConstIt secIt = firSource.begin();
        for (unsigned ind = 0; secIt != firSource.end(); ++secIt, ++ind) {
            peList.push_back(LeafPathElement(*secIt, firDir, ind));
        }
        return;
    }
    LPECollectionConstIt firIt = firSource.begin();
    LPECollectionConstIt secIt = secSource.begin();
    unsigned firInd = 0;
    unsigned secInd = 0;
    
    assert(peList.empty());
    bool tiesGoToFirst = (*secIt < *firIt ? false : true);
    bool secEmpty = false;
    bool firEmpty = false;
    while (peList.size() < gMaxNumLeafPaths && (!firEmpty || !secEmpty)) {
        const LeafPathElement *nLPE = 0L;
        TreeDirection d = LEFT_DIR;
        unsigned nextInd = 0;
        if (tiesGoToFirst) {
            // default to firIt in ties
            if ((!secEmpty) && (firEmpty || (*secIt < *firIt))) {
                nLPE = &(*secIt++);
                d = secDir;
                nextInd = secInd++;
                tiesGoToFirst = true;
                if (secIt == secSource.end()) {
                    secEmpty = true;
                }
            }
            else {
                assert(!firEmpty);
                nLPE = &(*firIt++);
                d = firDir;
                nextInd = firInd++;
                tiesGoToFirst = false;
                if (firIt == firSource.end()) {
                    firEmpty = true;
                }
            }
        }
        else {
            // default to secIt in ties
            if ((!firEmpty) && (secEmpty || (*firIt < *secIt))) {
                nLPE = &(*firIt++);
                d = firDir;
                nextInd = firInd++;
                tiesGoToFirst = false;
                if (firIt == firSource.end()) {
                    firEmpty = true;
                }
            }
            else {
                assert(!secEmpty);
                nLPE = &(*secIt++);
                d = secDir;
                nextInd = secInd++;
                tiesGoToFirst = true;
                if (secIt == secSource.end()) {
                    secEmpty = true;
                }
            }
        }

        assert(nLPE);
        peList.push_back(LeafPathElement(*nLPE, d, nextInd));
    }
    
    //std::cerr << "mergePathElementLists size = " << peList.size() << std::endl;

}


inline int CentroidScore(const NxsSimpleNode &nd,  long numActive) {
    NdBlob * b = (NdBlob *)nd.scratch;
    assert(b);
    const int above = b->numActiveLeavesAboveEdge;
    const int below = numActive - above;
    
    return (above > below ? above - below : below - above);
}
inline void mergeShortLeavesLists(const NxsSimpleNode &nd, TreeSweepDirection dir, const NxsSimpleNode &o, const NxsSimpleNode &o2) {
    if (gIntercalate) {
        assert(false);
        throw 5;
    }
    else {
        EdgeDecompInfo * edgeInfo = ((NdBlob*)nd.scratch)->activeEdgeInfo;
        NdBlob* oBlob = (NdBlob*)(o.scratch);
        NdBlob* o2Blob = (NdBlob*)(o2.scratch);
        assert(edgeInfo);
        if (dir == LEFT_OR_RIGHT) {
            //std::cerr << "mergeShortLeavesLists nd = " << (long) &nd << " ndblob = " << (long) nd.scratch << " ndBlob->activeEdgeInfo = " << (long)edgeInfo << " edgeInfo->closestLeavesAbove.size() = " << edgeInfo->closestLeavesAbove.size() << std::endl;
            //std::cerr << "mergeShortLeavesLists o = " << (long) &o << " oblob = " << (long) oBlob << " oBlob->activeEdgeInfo = " << (long)oBlob->activeEdgeInfo << " oBlob->activeEdgeInfo->closestLeavesAbove.size() = " << oBlob->activeEdgeInfo->closestLeavesAbove.size() << std::endl;
            //std::cerr << "mergeShortLeavesLists o22 = " << (long) &o2 << " o2blob = " << (long) o2Blob << " o2Blob->activeEdgeInfo = " << (long)oBlob->activeEdgeInfo << " oBlob->activeEdgeInfo->closestLeavesAbove.size() = " << o2Blob->activeEdgeInfo->closestLeavesAbove.size() << std::endl;
            LPECollection &peList = edgeInfo->closestLeavesAbove;
            assert(peList.empty());
            mergePathElementLists(peList, LEFT_DIR, oBlob->activeEdgeInfo->closestLeavesAbove,
                                          RIGHT_DIR, o2Blob->activeEdgeInfo->closestLeavesAbove);
        }
        else {
            const NxsSimpleNode * belowNd = &o2;
            EdgeDecompInfo * parInfo = o2Blob->activeEdgeInfo;
            EdgeDecompInfo * oInfo = oBlob->activeEdgeInfo;
            bool currIsLeft = (&nd == LeftChild(*belowNd));
            const LPECollection & parSource = parInfo->closestLeavesBelow;
            LPECollection &peList = edgeInfo->closestLeavesBelow;
            TreeDirection od = (dir == LEFT_OR_BELOW ? LEFT_DIR : RIGHT_DIR);
            mergePathElementLists(peList, od, oInfo->closestLeavesAbove, BELOW_DIR, parSource);
        }
    }
}

void nonRecursiveFlagActivePathUp(const NxsSimpleNode * currAnc, 
                                  const LPECollection * newActiveLeaves,
                                  const NxsSimpleNode * & nextNd,
                                  const LPECollection * & nextLPEC,
                                  const NxsSimpleNode * & futureNd,
                                  const LPECollection * & futureLPEC) {
    if (gDebugging) {
        std::cerr << "nonRecursiveFlagActivePathUp(currAnc=" << (long)currAnc << ", newActiveLeaves=" << (long) newActiveLeaves << ')' << std::endl;
        std::cerr << "LeftChild(*currAnc)=" << (long)currAnc->GetFirstChild() << std::endl;
        if (currAnc->GetFirstChild())
            std::cerr << "RightChild(*currAnc)=" << (long)RightChild(*currAnc)  << std::endl;
    }
    assert(currAnc);
    NdBlob * currBlob = (NdBlob*)(currAnc->scratch);
    assert(currBlob);
    
    currBlob->activeLeafDirStack.push(currBlob->activeLeafDir);
    currBlob->edgeInfoStack.push(currBlob->activeEdgeInfo);
    currBlob->activeEdgeInfo = 0L;
    currBlob->activeLeafDir = THIS_NODE_DIR_BIT;
    
    currBlob->lpeScratch1.clear();
    currBlob->lpeScratch2.clear();
    nextNd = 0L;
    nextLPEC = 0L;
    futureNd = 0L;
    futureLPEC = 0L;
    if (currAnc->IsTip()) {
        return;
    }
    if (newActiveLeaves == 0L) {
        LPECollectionConstIt nalIt = currBlob->fullEdgeInfo.closestLeavesAbove.begin();
        currBlob->activeEdgeInfo = &(currBlob->fullEdgeInfo);
        currBlob->numActiveLeavesAboveEdge = currBlob->fullEdgeInfo.closestLeavesAbove.size();
        for (; nalIt != currBlob->fullEdgeInfo.closestLeavesAbove.end(); ++nalIt) {
            const LeafPathElement & el = *nalIt;
            if (el.dirToNext == LEFT_DIR) {
                currBlob->lpeScratch1.push_back(el);
            }
            else {
                assert(el.dirToNext == RIGHT_DIR);
                currBlob->lpeScratch2.push_back(el);
            }
        }
    
    }
    else {
        LPECollectionConstIt nalIt = newActiveLeaves->begin();
        currBlob->activeEdgeInfo = new EdgeDecompInfo();
        currBlob->numActiveLeavesAboveEdge = newActiveLeaves->size();
        for (; nalIt != newActiveLeaves->end(); ++nalIt) {
            int indInCurr = nalIt->indexInNext;
            const LeafPathElement & el = currBlob->fullEdgeInfo.closestLeavesAbove.at(indInCurr);
            if (gDebugging) {
                std::cerr << "nalIt->indexInNext = " << indInCurr << " el.leaf = " << (long) el.leaf << " nalIt->leaf = " << (long) nalIt->leaf << std::endl;
            }
            assert(el.leaf == nalIt->leaf);
            if (el.dirToNext == LEFT_DIR) {
                currBlob->lpeScratch1.push_back(el);
            }
            else {
                assert(el.dirToNext == RIGHT_DIR);
                currBlob->lpeScratch2.push_back(el);
            }            
        }
    }
    if (currBlob->lpeScratch1.empty()) {
        if (currBlob->lpeScratch2.empty()) {
            return;
        }
        currBlob->activeLeafDir = RIGHT_DIR_BIT;
        futureNd = RightChild(*currAnc);
        futureLPEC = &(currBlob->lpeScratch2);
    }
    else if (currBlob->lpeScratch2.empty()) {
        currBlob->activeLeafDir = LEFT_DIR_BIT;
        nextNd = LeftChild(*currAnc);
        nextLPEC = &(currBlob->lpeScratch1);
        }
    else {
        currBlob->activeLeafDir = LEFT_OR_RIGHT;
        nextNd = LeftChild(*currAnc);
        nextLPEC = &(currBlob->lpeScratch1);
        futureNd = RightChild(*currAnc);
        futureLPEC = &(currBlob->lpeScratch2);
    }
    
}

// Move up the tree..
// Pushing current activeEdgeInfo and activeLeafDir to their stacks and update based
//  on the `newActiveLeaves`  If `newActiveLeaves` is 0L, then the entire subtree
//  should be regarded as active
void flagActivePathUp(const NxsSimpleNode *currAnc, const LPECollection *newActiveLeaves) {
    if (gDebugging) {
        std::cerr << "flagActivePathUp(currAnc=" << (long) currAnc << ", newActiveLeaves=" <<(long)newActiveLeaves << ")" << std::endl;
        if (newActiveLeaves) {
            std::cerr << " leaf paths:\n";
            for (LPECollectionConstIt i = newActiveLeaves->begin(); i != newActiveLeaves->end(); ++i) {
                std::cerr << "  " << *i << '\n';
            }
            std::cerr << std::endl;
        }
    }
    std::stack<const NxsSimpleNode *> ndStack;
    std::stack<const LPECollection *> lpecStack;
    ndStack.push(currAnc);
    lpecStack.push(newActiveLeaves);
    while (!ndStack.empty()) {
        const NxsSimpleNode * currNd = ndStack.top();
        ndStack.pop();
        assert(!lpecStack.empty());
        const LPECollection * currLPEC = lpecStack.top();
        lpecStack.pop();
        
        const NxsSimpleNode * nextNd = 0L;
        const NxsSimpleNode * futureNd = 0L;
        const LPECollection * nextLPEC = 0L;
        const LPECollection * futureLPEC = 0L;
        nonRecursiveFlagActivePathUp(currNd, currLPEC, nextNd, nextLPEC, futureNd, futureLPEC);
        if (futureNd != 0) {
            assert(futureLPEC != 0);
            ndStack.push(futureNd);
            lpecStack.push(futureLPEC);
        }
        if (nextNd != 0) {
            assert(nextLPEC != 0);
            ndStack.push(nextNd);
            lpecStack.push(nextLPEC);
        }
    }
}
void flagActivePathDown(const NxsSimpleNode *currAnc, const LPECollection *newActiveLeaves) {
    std::stack<const NxsSimpleNode *> ndStack;
    std::stack<const LPECollection *> collStack;

}

////////////////////////////////////////////////////////////////////////////////
// leftSubtreeRoot=====>v   v <=rightSubtreeRoot
//                       \ /
//   sibSubtreeRoot=> v   o <=centroidChild
//                     \ / 
//                par=> o   v
//                       \ /
//                    v   o <= gpSubtreeRoot (grandparent)
//                     \ /
//                      o
//                     /
//##############################################################################
// Or (if gpSubtreeRoot would be the root of the tree):
//##############################################################################
//        leftSubtreeRoot=====>v   v <=rightSubtreeRoot
//                              \ /
//          sibSubtreeRoot=> v   o <=centroidChild
//                            \ / 
// (leaf0, gpSubtreeRoot)=>o   o
//                          \ /
//                           o <= (iniRoot of the tree)
///////////////////////////////////////////////////////////////////////////
void DecomposeAroundCentroidChild(std::vector<const NxsSimpleNode *> &preorderTraversal,
                                  const NxsSimpleNode *centroidChild,
                                  long numActiveLeaves) {
    ////////////////////////////////////////////////////////////////////////////
    // Step one is to get identify leftSubtreeRoot, rightSubtreeRoot, sibSubtreeRoot
    //  and gpSubtreeRoot (see diagrams above).
    ////////////////////////////////////////////////////////////////////////////
    assert(centroidChild);
    assert(!centroidChild->IsTip());
    const NxsSimpleNode * leftSubtreeRoot = LeftChild(*centroidChild);
    assert(leftSubtreeRoot);
    const NxsSimpleNode * rightSubtreeRoot = RightChild(*centroidChild);
    assert(rightSubtreeRoot);
    const NxsSimpleNode * par = Parent(*centroidChild);
    assert(par);
    assert(!par->IsTip());
    const NxsSimpleNode * parlc = LeftChild(*par);
    const NxsSimpleNode * sibSubtreeRoot = 0L;
    if (parlc == centroidChild) {
        sibSubtreeRoot = parlc->GetNextSib();
        assert(sibSubtreeRoot);
        assert(sibSubtreeRoot->GetNextSib() == 0L);
    }
    else {
        assert(parlc->GetNextSib() == centroidChild);
        assert(centroidChild->GetNextSib() == 0L);
        sibSubtreeRoot = parlc;
    }
    assert(sibSubtreeRoot);
    
    const NxsSimpleNode * gpSubtreeRoot = Parent(*par);
    assert(gpSubtreeRoot);
    if (IsRootChild(*par)) {
        gpSubtreeRoot = LeftChild(*Parent(*par));
        assert(gpSubtreeRoot->GetNextSib() == par);
        assert(par->GetNextSib() == 0L);
    }
    else {
        assert(Parent(*gpSubtreeRoot));
    }
    
    NdBlob * centroidBlob = (NdBlob *)(centroidChild->scratch);
    NdBlob * gpBlob = (NdBlob *)(par->scratch); // Note that parBlob is actually from the edge that connects par to gpSubtreeRoot
    NdBlob * sibBlob = (NdBlob *)(sibSubtreeRoot->scratch); // Note that parBlob is actually from the edge that connects par to gpSubtreeRoot
    NdBlob * leftBlob = (NdBlob *)(leftSubtreeRoot->scratch);
    NdBlob * rightBlob = (NdBlob *)(rightSubtreeRoot->scratch);
    ////////////////////////////////////////////////////////////////////////////
    // Step two is to identify the common leaf set for all four decompositions.
    ////////////////////////////////////////////////////////////////////////////
    long numActLeavesAboveCentroid = centroidBlob->numActiveLeavesAboveEdge;
    long numActLeavesBelowCentroid = numActiveLeaves - numActLeavesAboveCentroid;
    long numAboveChosen, numBelowChosen, numLeftChosen, numRightChosen, numSibChosen, numGPChosen;
    if (numActLeavesBelowCentroid < gFloorHalfIntersectionSize) {
        numBelowChosen = numActLeavesBelowCentroid;
        assert(numActLeavesAboveCentroid >= gLeafSetIntersectionSize - numBelowChosen);
        numAboveChosen = gLeafSetIntersectionSize - numBelowChosen;
    }
    else {
        if (numActLeavesAboveCentroid < gFloorHalfIntersectionSize) {
            numAboveChosen = numActLeavesAboveCentroid;
            assert(numActLeavesBelowCentroid >= gLeafSetIntersectionSize - numAboveChosen);
        }
        else {
            numAboveChosen = gFloorHalfIntersectionSize; //\TODO should probably sort the PLE's and use their order to decide which direction gets the rounding error...
        }
        numBelowChosen = gLeafSetIntersectionSize - numAboveChosen;
    }

    long floorHalfAbove = numAboveChosen/2;
    if (leftBlob->numActiveLeavesAboveEdge < floorHalfAbove) {
        numLeftChosen = leftBlob->numActiveLeavesAboveEdge;
        assert(rightBlob->numActiveLeavesAboveEdge >= floorHalfAbove - numLeftChosen);
        numRightChosen = numAboveChosen - numLeftChosen;
    }
    else {
        if (rightBlob->numActiveLeavesAboveEdge < floorHalfAbove) {
            numRightChosen = rightBlob->numActiveLeavesAboveEdge;
            assert(leftBlob->numActiveLeavesAboveEdge >= numAboveChosen - numRightChosen);
            
        }
        else {
            numRightChosen = floorHalfAbove; //\TODO should probably sort the PLE's and use their order to decide which direction gets the rounding error...
        }
        numLeftChosen = numAboveChosen - numRightChosen;
    }
    

    long floorHalfBelow = numBelowChosen/2;
    if (sibBlob->numActiveLeavesAboveEdge < floorHalfBelow) {
        numSibChosen = sibBlob->numActiveLeavesAboveEdge;
        assert((numActiveLeaves  - gpBlob->numActiveLeavesAboveEdge) >= numBelowChosen - numSibChosen);
        numGPChosen = numBelowChosen - numSibChosen;
    }
    else {
        if ((numActiveLeaves  - gpBlob->numActiveLeavesAboveEdge) < floorHalfBelow) {
            numGPChosen = (numActiveLeaves  - gpBlob->numActiveLeavesAboveEdge);
            assert(sibBlob->numActiveLeavesAboveEdge >= numBelowChosen - numGPChosen);
        }
        else {
            numGPChosen = floorHalfBelow; //\TODO should probably sort the PLE's and use their order to decide which direction gets the rounding error...
        }
        numSibChosen = numBelowChosen - numGPChosen;
    }

    assert(numAboveChosen + numBelowChosen == gLeafSetIntersectionSize);
    assert(numLeftChosen + numRightChosen == numAboveChosen);
    assert(numSibChosen + numGPChosen == numBelowChosen);
    assert(numLeftChosen > 0);
    assert(numRightChosen > 0);
    assert(numSibChosen > 0);
    assert(numGPChosen > 0);
    
    
    const EdgeDecompInfo * leftEDI = leftBlob->activeEdgeInfo;
    assert(leftEDI);
    const LPECollection & leftFullLPEC = leftEDI->closestLeavesAbove;
    assert(leftFullLPEC.size() >= numLeftChosen);
    LPECollection::const_iterator lfLPECIt = leftFullLPEC.begin();
    LPECollection leftCommonLeafSet;
    leftCommonLeafSet.reserve(numLeftChosen);
    for (unsigned i=0; i < numLeftChosen; ++i, ++lfLPECIt) {
        leftCommonLeafSet.push_back(LeafPathElement(*lfLPECIt, LEFT_DIR, i));
    }
    assert(leftCommonLeafSet.size() == numLeftChosen); 
    
    const EdgeDecompInfo * rightEDI = rightBlob->activeEdgeInfo;
    assert(rightEDI);
    const LPECollection & rightFullLPEC = rightEDI->closestLeavesAbove;
    assert(rightFullLPEC.size() >= numRightChosen);
    LPECollection::const_iterator rfLPECIt = rightFullLPEC.begin();
    LPECollection rightCommonLeafSet;
    rightCommonLeafSet.reserve(numRightChosen);
    for (unsigned i=0; i < numRightChosen; ++i, ++rfLPECIt) {
        rightCommonLeafSet.push_back(LeafPathElement(*rfLPECIt, RIGHT_DIR, i));
    }
    assert(rightCommonLeafSet.size() == numRightChosen); 

    const EdgeDecompInfo * sibEDI = sibBlob->activeEdgeInfo;
    assert(sibEDI);
    const LPECollection & sibFullLPEC = sibEDI->closestLeavesAbove;
    assert(sibFullLPEC.size() >= numSibChosen);
    LPECollection::const_iterator sibfLPECIt = sibFullLPEC.begin();
    LPECollection sibCommonLeafSet;
    sibCommonLeafSet.reserve(numSibChosen);
    for (unsigned i=0; i < numSibChosen; ++i, ++sibfLPECIt) {
        sibCommonLeafSet.push_back(LeafPathElement(*sibfLPECIt, (sibBlob->isParentsLeftChild ? LEFT_DIR : RIGHT_DIR), i));
    }
    assert(sibCommonLeafSet.size() == numSibChosen); 

    const EdgeDecompInfo * gpEDI = gpBlob->activeEdgeInfo;
    assert(gpEDI);
    const LPECollection & gpFullLPEC = gpEDI->closestLeavesBelow;
    assert(gpFullLPEC.size() >= numGPChosen);
    LPECollection::const_iterator gpfLPECIt = gpFullLPEC.begin();
    LPECollection gpCommonLeafSet;
    gpCommonLeafSet.reserve(numGPChosen);
    for (unsigned i=0; i < numGPChosen; ++i, ++gpfLPECIt) {
        gpCommonLeafSet.push_back(LeafPathElement(*gpfLPECIt, BELOW_DIR, i));
    }
    assert(gpCommonLeafSet.size() == numGPChosen); 
    
    ////////////////////////////////////////////////////////////////////////////
    // Step 3 - 6 activate each subproblem in turn and recurse...
    //
    ////////////////////////////////////////////////////////////////////////////
    flagActivePathUp(leftSubtreeRoot, 0L);
    flagActivePathUp(rightSubtreeRoot, &rightCommonLeafSet);
    flagActivePathUp(sibSubtreeRoot, &sibCommonLeafSet);
    flagActivePathDown(par, &gpCommonLeafSet);
    long sizeOfLeftSubproblem = numBelowChosen + numRightChosen + leftBlob->numActiveLeavesAboveEdge;


    
    
}


bool treeReadCallback(NxsFullTreeDescription &ftd, void *x, NxsTreesBlock *treesBlock) {
    assert(treesBlock != 0L);
    assert(gOutputStream != 0L);
    if (!gTaxaBlockWritten) {
        NxsTaxaBlock * taxa = (NxsTaxaBlock *)treesBlock->GetTaxaBlockPtr();
        if (gPrevTaxaBlock != 0L) {
            if (gPrevTaxaBlock != taxa) {
                gErrorMessage << "A second taxon block was encountered. Multiple taxa blocks are not supported." ;
                throw NxsException(gErrorMessage);
            }
        }
        else 
            gPrevTaxaBlock = taxa;
        if (gOutputStream) {
            *gOutputStream << "#NEXUS\n";
            if (taxa != 0L) {
                taxa->WriteAsNexus(*gOutputStream);                
            }
            *gOutputStream << "BEGIN TREES;\n";
            treesBlock->WriteTranslateCommand(*gOutputStream);
        }
        gTaxaBlockWritten = true;
        gTaxaBlockNumLeaves = (long) taxa->GetNTax();
        if (gTaxaBlockNumLeaves < gLeafSetIntersectionSize) {
            gErrorMessage << "LeafSetIntersectionSize = " << gLeafSetIntersectionSize << ", but there are only " << gTaxaBlockNumLeaves << " leaves in the input file." ;
            throw NxsException(gErrorMessage);
        }
    }
    if (ftd.HasPolytomies()) {
        gErrorMessage << "Tree # " << gCurrTreeIndex + 1 << " contains a polytomy!";
        throw NxsException(gErrorMessage);
    }
    if (gUseEdgeLengths && !ftd.AllEdgesHaveLengths()) {
        gErrorMessage << "The use EDGE_LENGTHS option is in effect, but tree " << gCurrTreeIndex + 1 << " contains edges without lengths!";
        throw NxsException(gErrorMessage);
    }
    if (ftd.HasDegreeTwoNodes()) {
        gErrorMessage << "Tree " << gCurrTreeIndex + 1 << " has nodes of outdegree 1!";
        throw NxsException(gErrorMessage);
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // To simplify processing of the tree we make the following modifications 
    // of the NxsSimpleTree:
    //      1. make leaf 0 the left child of the root
    //      2. make the root have two children (introduce bogusNode as the root
    //          if the original root (iniRoot) has three children).
    //      3. call the right child of the root the "virtual root" or (virtRoot)
    //      4. prepare a vector of a nodes that is the preorder traversal of the
    //          subtree rooted at virtRoot:
    //
    //                 __ __
    //                 \/ \/
    //                  \ /
    //       leaf0=> o   o <=virtRoot
    //                \ /
    //       iniRoot=> o
    ////////
    NxsSimpleTree t(ftd, 0, 0.0);
    NxsSimpleNode * leaf0 = t.RerootAt(0);
    NxsSimpleNode * iniRoot = const_cast<NxsSimpleNode *>(t.GetRootConst());
    NxsSimpleNode * rootLeft = const_cast<NxsSimpleNode *>(LeftChild(*iniRoot));
    if (rootLeft != leaf0) {
        if (rootLeft->GetNextSib() == leaf0) {
            rootLeft->LowLevelSetNextSib(leaf0->GetNextSib());   
        }
        else { 
            if ((rootLeft->GetNextSib() == 0) || (rootLeft->GetNextSib()->GetNextSib() != leaf0) || (leaf0->GetNextSib() != 0)) {
                gErrorMessage << "Rerooting at leaf 0 did not succeed!";
                throw NxsException(gErrorMessage);
            }
            NxsSimpleNode * rootMiddle = rootLeft->GetNextSib()->GetNextSib();
            rootMiddle->LowLevelSetNextSib(0L);
        }
        leaf0->LowLevelSetNextSib(rootLeft);
        iniRoot->LowLevelSetFirstChild(leaf0);
    }
    NxsSimpleNode * toDelNd = 0L;
    std::vector<NxsSimpleNode *> rootChildren = iniRoot->GetChildren();
    gRoot = iniRoot;
    if (rootChildren.size() != 2) {
        if (rootChildren.size() < 2) {
            gErrorMessage << "The root of tree " << gCurrTreeIndex + 1 << " has fewer than 2 children!";
            throw NxsException(gErrorMessage);
        }
        NxsSimpleNode * rlc = rootChildren.at(0);
        iniRoot->RemoveChild(rlc);

        bogusNode.AddChild(rlc);
        NxsSimpleEdge & e1 = rlc->GetMutableEdgeToParentRef();
        rlc->LowLevelSetNextSib(0L);
        e1.SetParent(&bogusNode);

        bogusNode.AddChild(iniRoot);
        NxsSimpleEdge & e2 = iniRoot->GetMutableEdgeToParentRef();
        e2.SetParent(&bogusNode);

        gRoot = & bogusNode;
    }
    
    std::vector<const NxsSimpleNode *> preorderTraversal;
    NxsSimpleNode * virtRoot = rootLeft->GetNextSib();
    assert(virtRoot);
    assert(virtRoot->GetNextSib() == 0L);
    virtRoot->AddSelfAndDesToPreorder(preorderTraversal);

    gNumNodes = (long) preorderTraversal.size() + 2;
    if ((gNumNodes % 2) == 0) {
        gErrorMessage << "Tree " << gCurrTreeIndex + 1 << " has an even number of nodes!";
        throw NxsException(gErrorMessage);
    }
    gNumLeaves = (long)((gNumNodes + 1 ) / 2);
    if (gNumLeaves < gLeafSetIntersectionSize) {
        gErrorMessage << "LeafSetIntersectionSize = " << gLeafSetIntersectionSize << ", but there are only " << gNumLeaves << " leaves in tree " << gCurrTreeIndex + 1  << '.';
        throw NxsException(gErrorMessage);
    }
    
    std::vector<NdBlob *> gAllocedBlobs;
    gAllocedBlobs.resize(gNumNodes);
    std::vector<NxsSimpleNode *> bogusChildren = bogusNode.GetChildren();
    
    try {
        const NxsSimpleNode * centroidChild = 0L;
        int minCentroidScore = gNumLeaves;
        
        if (gNumLeaves <= gMaxSubProblemSize) {
            if (gOutputStream) {
                *gOutputStream << "Tree tree" << gCurrTreeIndex << " = " << ftd.GetNewick() << ";" << std::endl;
            }
        } 
        else {
            
            NdBlob * rlb = new NdBlob(true);
            rlb->isParentsLeftChild = true;
            rootLeft->scratch = (void *) rlb;
            rlb->numLeavesAboveEdge = 1;
            rlb->numActiveLeavesAboveEdge = 1;
            EdgeDecompInfo * rlEdgeInfo = rlb->activeEdgeInfo;
            rlEdgeInfo->closestLeavesAbove.push_back(LeafPathElement(rootLeft, THIS_NODE_DIR));
            rlb->activeLeafDir = ALL_DIR_BITS;
            
            std::vector<const NxsSimpleNode *>::const_reverse_iterator postIt = preorderTraversal.rbegin();
            
            for (; postIt != preorderTraversal.rend(); ++postIt) {
                const NxsSimpleNode * nd = *postIt;
                NdBlob * nb = new NdBlob(LeftChild(*Parent(*nd)) == nd);
                gAllocedBlobs.push_back(nb);
                nd->scratch = (void *) nb;
                nb->activeLeafDir = ALL_DIR_BITS;

                if (nd->IsTip()) {
                    nb->numLeavesAboveEdge = 1;
                    nb->numActiveLeavesAboveEdge = 1;
                    EdgeDecompInfo * edgeInfo = nb->activeEdgeInfo;
                    edgeInfo->closestLeavesAbove.push_back(LeafPathElement(nd, THIS_NODE_DIR));
                }
                else {
                    nb->numLeavesAboveEdge = LeftBlob(*nd)->numLeavesAboveEdge + RightBlob(*nd)->numLeavesAboveEdge;
                    nb->numActiveLeavesAboveEdge = nb->numLeavesAboveEdge;
                    int currSC = CentroidScore(*nd, gNumLeaves);
                    if (currSC < minCentroidScore) {
                        centroidChild = nd;
                        minCentroidScore = currSC;
                        if (gDebugging) {
                            std::cerr << "New Min Centroid Score = " << minCentroidScore << std::endl;
                        }
                    }
                    mergeShortLeavesLists(*nd, LEFT_OR_RIGHT, *LeftChild(*nd), *RightChild(*nd));
                    
                }
            }
            NdBlob * rb = new NdBlob(false);
            rb->isParentsLeftChild = false;
            gRoot->scratch = (void *) rb;
            rb->numLeavesAboveEdge = LeftBlob(*gRoot)->numLeavesAboveEdge + RightBlob(*gRoot)->numLeavesAboveEdge;
            rb->numActiveLeavesAboveEdge = rb->numLeavesAboveEdge;
            rb->activeLeafDir = ALL_DIR_BITS;

            mergeShortLeavesLists(*gRoot, LEFT_OR_RIGHT, *LeftChild(*gRoot), *RightChild(*gRoot));
                    
            if (((NdBlob *)gRoot->scratch)->numLeavesAboveEdge != gNumLeaves) {
                gErrorMessage << "Expected " << gNumLeaves << " leaves at the root, but only found " << ((NdBlob *)gRoot->scratch)->numLeavesAboveEdge;
                throw NxsException(gErrorMessage);
            }
            
            std::vector<const NxsSimpleNode *>::const_iterator preIt = preorderTraversal.begin();
            preIt++; //skip the root
            for (; preIt != preorderTraversal.end(); ++preIt) {
                const NxsSimpleNode * nd = *preIt;
                assert(nd != 0L);
                if (!nd->IsTip()) {
                    NdBlob* ndBlob = (NdBlob*)(nd->scratch);
                    EdgeDecompInfo * edgeInfo = ndBlob->activeEdgeInfo;
                    assert(edgeInfo);
                    if (IsRootChild(*nd)) {
                        const NxsSimpleNode * leafChild = gRoot->GetFirstChild();
                        assert(leafChild != nd);
                        assert(leafChild->GetNextSib() == nd);
                        assert(leafChild->IsTip());
                        EdgeDecompInfo * lcEdgeInfo = ((NdBlob*)leafChild->scratch)->activeEdgeInfo;
                        assert(lcEdgeInfo);
                        assert(lcEdgeInfo->closestLeavesAbove.size() == 1);
                        assert(edgeInfo->closestLeavesBelow.empty());
                        edgeInfo->closestLeavesBelow.push_back(LeafPathElement(lcEdgeInfo->closestLeavesAbove[0], BELOW_DIR, 0));
                    }
                    else {
                        if (ndBlob->isParentsLeftChild)
                            mergeShortLeavesLists(*nd, RIGHT_OR_BELOW, *RightChild(*nd), *Parent(*nd));
                        else
                            mergeShortLeavesLists(*nd, LEFT_OR_BELOW, *LeftChild(*nd), *Parent(*nd));
                    }
                }
                if (gDebugging) {
                    NdBlob * nb = (NdBlob *) nd->scratch;
                    std::cerr << "nd = " << (long) nd << " blob=" << (long) nb << " blob->numLeavesAboveEdge = " << nb->numLeavesAboveEdge <<  " blob->activeEdgeInfo = " << (long) nb->activeEdgeInfo <<  " blob->activeEdgeInfo->closestLeavesAbove.size() = " << nb->activeEdgeInfo->closestLeavesAbove.size() << '\n';
                    if (nd->IsTip()) {
                        std::cerr << "nd is tip, so below info not calculated" << std::endl;
                    }
                    else {
                        std::cerr << "nd = " << (long) nd << " blob=" << (long) nb << " blob->numLeavesBelowEdge = " << (gNumLeaves - nb->numLeavesAboveEdge) <<  " blob->activeEdgeInfo = " << (long) nb->activeEdgeInfo <<  " blob->activeEdgeInfo->closestLeavesBelow.size() = " << nb->activeEdgeInfo->closestLeavesBelow.size() << '\n';
                    }
                }
            }            
        }        
        
        assert(centroidChild != 0);
        
        DecomposeAroundCentroidChild(preorderTraversal, centroidChild, gNumLeaves);
    }
    catch (...) {
         for (std::vector<NxsSimpleNode *>::iterator cIt = bogusChildren.begin(); cIt != bogusChildren.end(); ++cIt)
            bogusNode.RemoveChild(*cIt);
        for (std::vector<NdBlob *>::iterator ab = gAllocedBlobs.begin(); ab != gAllocedBlobs.end(); ++ab)
            delete *ab;
        throw;
        }

    for (std::vector<NxsSimpleNode *>::iterator cIt = bogusChildren.begin(); cIt != bogusChildren.end(); ++cIt)
        bogusNode.RemoveChild(*cIt);
    for (std::vector<NdBlob *>::iterator ab = gAllocedBlobs.begin(); ab != gAllocedBlobs.end(); ++ab)
        delete *ab;
    gCurrTreeIndex++;
    return false;
}


int readInput(std::istream &inp, MultiFormatReader::DataFormatType fmt, const std::string & filename) {
    assert(gNexusReader == 0L);
    int blockFlag = PublicNexusReader::NEXUS_TAXA_BLOCK_BIT | PublicNexusReader::NEXUS_TREES_BLOCK_BIT;
	gNexusReader = new MultiFormatReader(blockFlag, NxsReader::WARNINGS_TO_STDERR);
	if (!gVerbose)
		gNexusReader->SetWarningOutputLevel(NxsReader::SKIPPING_CONTENT_WARNING);
	NxsTreesBlock * treesB = gNexusReader->GetTreesBlockTemplate();
	assert(treesB);
	treesB->SetAllowImplicitNames(true);
    treesB->setValidationCallbacks(treeReadCallback, 0L);
    treesB->SetTreatAsRootedByDefault(false);
    try {
        gNexusReader->ReadStream(inp, fmt, filename.c_str());
    }
    catch (NxsException & x) {
		std::cerr << "Error:\n " << x.msg << std::endl;
		std::cerr << "File: \"" << filename << '\"' << std::endl;
		if (x.line > 0 || x.pos > 0)
			std::cerr << "At line " << x.line << ", column (approximately) " << x.col << " (and file position "<< x.pos << ")" << std::endl;
		return 14;
    }
    catch (...) {
        std::cerr << gErrMsgPrefix << "Unknown exception generated when reading \"" << filename << '\"' << std::endl;
		return 15;
    }
    return 0;
}

void printHelp(std::ostream & out)
	{
	out << gVersionString << "\n";
	out << "Reads a file with a tree description and produces a decomposition of the tree\n";
	out << "into smaller trees.  Subject to the following constraints:\n";
	out << "    1. No tree will have more than MAX_SUBTREE_SIZE leaves\n";
	out << "    2. Trees will be written in hierarchichal groups (with their names\n";
	out << "           indicating their group membership)\n";
	out << "    3. Within each group, each tree will contain a common subset of leaves.\n";
	out << "           The size of this subset will be at least MIN_LEAF_INTERSECTION_SIZE.\n";
	out << "    4. All leaves in the original tree will be present in at least one output tree.\n";
	out << "\nThe input tree must be fully-resolved (binary).\n";
	out << "\nDecomposition is produced by recursively breaking the tree around centroid\n";
	out << "    edge and striving to include a group of MIN_LEAF_INTERSECTION_SIZE/4 leaves\n";
	out << "    from each induced subtree into the smaller tree (if some subtrees are smaller\n";
	out << "    than MIN_LEAF_INTERSECTION_SIZE/4 then more leaves will be included from their\n";
	out << "    \"sister\" subtree).\n";
	out << "\nAt each step, leaves closest to the decomposition edge are retained.  This distance.\n";
	out << "    can be topological or based on EDGE_LENGTHS (-e option below; in this case all edges\n";
	out << "    in the tree must have edge lengths).\n";
	out << "\nCommand-line flags:\n\n";
	out << "    -d debugging output.\n\n";
	out << "    -e use EDGE_LENGTHS in the shortest leaf calculations.\n\n";
	out << "    -f<format> specifies the input file format expected:\n";
	out << "            -fnexus     NEXUS (this is also the default)\n";
	out << "            -fphyliptree    Phylip tree file (simple newick string)\n";
	out << "            -frelaxedphyliptree Relaxed Phylip name restrictions\n";
	out << "    -h on the command line shows this help message.\n\n";
	out << "    -i<non-negative integer> specifies MIN_LEAF_INTERSECTION_SIZE\n";
	out << "    -m<non-negative integer> specifies MAX_SUBTREE_SIZE\n";
	out << "    -v verbose output to stdandard error.\n\n";
}


int main(int argc, char * argv[]) {
	NxsReader::setNCLCatchesSignals(true);
	MultiFormatReader::DataFormatType f(MultiFormatReader::NEXUS_FORMAT);
    std::string formatName("NEXUS");
    
	for (int i = 1; i < argc; ++i) {
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (slen < 2 || filepath[0] != '-')
			continue;
		if (filepath[1] == 'h') {
			printHelp(std::cout);
			return 0;
		}
		else if (filepath[1] == 'd') {
		    gDebugging = true;
			gVerbose = true;
		}
		else if (filepath[1] == 'v')
			gVerbose = true;
		else if (filepath[1] == 'e')
			gUseEdgeLengths = true;
        else if (filepath[1] == 'i') {
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gLeafSetIntersectionSize))) {
				std::cerr << gErrMsgPrefix << "Expecting an integer -i\n" << std::endl;
				return 3;
			}
		}
        else if (filepath[1] == 'm') {
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &gMaxSubProblemSize))) {
				std::cerr << gErrMsgPrefix << "Expecting an integer -m\n" << std::endl;
				return 4;
			}
		}
        else if (filepath[1] == 'f') {
			f = MultiFormatReader::UNSUPPORTED_FORMAT;
			if (slen > 2) {
				formatName.assign(filepath + 2, slen - 2);
				f =  MultiFormatReader::formatNameToCode(formatName);
				if ((f != MultiFormatReader::NEXUS_FORMAT)
				    && (f != MultiFormatReader::PHYLIP_TREE_FORMAT)
				    && (f != MultiFormatReader::RELAXED_PHYLIP_TREE_FORMAT)
				   ){
					std::cerr << gErrMsgPrefix << "Unsupported tree file format \"" << formatName << "\" after -f\n" << std::endl;
					return 11;
                }
            }
			if (f == MultiFormatReader::UNSUPPORTED_FORMAT) {
				std::cerr << gErrMsgPrefix << "Expecting a format after -f\n" << std::endl;
				return 12;
            }
        }
        else if (strcmp(filepath, "--version") == 0) {
                std::cout << gVersionString << std::endl;
                return 0;
			}
		else {
            std::cerr << gErrMsgPrefix << "Unexpected command-line option: \"" << filepath<< "\"" << std::endl;
			return 5;
		}
	}
    if (gLeafSetIntersectionSize < 4) {
        std::cerr << gErrMsgPrefix << "Minimum size of the leaf intersection set must 4" << std::endl;
        return 6;
    }
    if (gMaxSubProblemSize < gLeafSetIntersectionSize + 3) {
        std::cerr << "Max Sub problem size must be at least 3 larger that the minimum size of the leaf intersection set" << std::endl;
        return 7;
    }


	bool filefound = false;
	std::string filename("<standard input>");
	for (int i = 1; i < argc; ++i) {
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (slen < 1){
            std::cerr << gErrMsgPrefix << "Unexpected empty string as command-line option " << std::endl;
			return 8;
		}
		if (filepath[0] == '-')
			continue;
        if (filefound) {
            std::cerr << gErrMsgPrefix << "Expecting only one file to read!" << std::endl;
			return 9;
        }
        filename.assign(filepath);
        filefound = true;
    }


	std::istream * inpStream = &(std::cin);
	std::ifstream inFileStream;


	if (filefound) {
	    inFileStream.open(filename.c_str());
	    if (!inFileStream.good()) {
            std::cerr << gErrMsgPrefix << "Could not open the file \""<< filename << '\"' << std::endl;
			return 10;
	    }
	    inpStream = &inFileStream;
	}

	if (gVerbose) {
	    // write invocation to std error, to make it easy to rerun...
	    std::cerr << argv[0] << " -v";
	    if (gUseEdgeLengths) {
	        std::cerr << " -e";
	    }
	    std::cerr << " -i" << gLeafSetIntersectionSize;
	    std::cerr << " -m" << gMaxSubProblemSize;
	    std::cerr << " -f" << formatName;
	    if (filefound) {
	        std::cerr << " \"" << filename  << "\"";
	    }
	    std::cerr << std::endl;
	}
	if (gVerbose && !filefound) {
	    std::cerr << gErrMsgPrefix << "Reading from stdin..." << std::endl;
	}

    gMaxNumLeafPaths = gLeafSetIntersectionSize; //@TODO this could be set to something smaller. gLeafSetIntersectionSize / 2, perhaps?

    gFloorHalfIntersectionSize = gLeafSetIntersectionSize/2;

    int rc =  readInput(*inpStream, f, filename);

    if (rc == 0 && gOutputStream != 0L) {
        *gOutputStream << "END;\n";
    }
    return rc;
}

