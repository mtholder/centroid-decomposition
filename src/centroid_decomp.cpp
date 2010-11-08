#include "ncl/nxsmultiformat.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <list>
#include <algorithm>

#include "centroid_decomp.hpp"


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
bool gEmitIntermediates = false;

NxsSimpleNode * rootLeft ; // leaf 0, which will be the left child of the root



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

inline void NdBlob::popBlob() {
    if (this->activeEdgeInfo != &(this->fullEdgeInfo)) {
        delete this->activeEdgeInfo;
    }
    this->activeEdgeInfo = this->edgeInfoStack.top();
    this->edgeInfoStack.pop();
    this->activeLeafDir  = this->activeLeafDirStack.top();
    this->activeLeafDirStack.pop();
    this->focalEdgeDir = this->focalEdgeDirStack.top();
    this->focalEdgeDirStack.pop();
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

inline std::ostream & operator<<(std::ostream & o, const TreeSweepDirection &d) {
    if (d == LEFT_DIR_BIT)
        o << "LEFT_DIR_BIT";
    else if (d == RIGHT_DIR_BIT)
        o << "RIGHT_DIR_BIT";
    else if (d == LEFT_OR_RIGHT)
        o << "LEFT_OR_RIGHT";
    else if (d == BELOW_DIR_BIT)
        o << "BELOW_DIR_BIT";
    else if (d == LEFT_OR_BELOW)
        o << "LEFT_OR_BELOW";
    else if (d == RIGHT_OR_BELOW)
        o << "RIGHT_OR_BELOW";
    else if (d == ALL_DIR_BITS)
        o << "ALL_DIR_BITS";
    else if (d == THIS_NODE_DIR_BIT)
        o << "THIS_NODE_DIR_BIT";
    else {
        assert(false);
        o << "Unknown direction";
        }
    return o;
}

inline std::ostream & operator<<(std::ostream & o, const LeafPathElement &lpe) {
    o << "LeafPathElement(leaf= "<< (long)lpe.leaf << ", dirToNext= " << lpe.dirToNext << ", indexInNext= " << lpe.indexInNext << ")";
    return o;
}


inline void writeEdgeDecompInfo(std::ostream & o, const EdgeDecompInfo & edi) {
    o << "edi.CLA = ";
    writeLeafSet(o, edi.closestLeavesAbove);
    o << "; edi.CLB = ";
    writeLeafSet(o, edi.closestLeavesBelow);
};

inline void debugBlobPrint(std::ostream & o, const NdBlob * b) {
    o << "b#LvsAbove = " << b->GetNumLeavesAbove();
    o << " b#ActLvsAbove = "<< b->GetNumActiveLeavesAbove();
    o << " b.actDir = "<< b->GetActiveLeafDir();
    o << " b.focalDir = "<< b->GetFocalEdgeDir();
    o << " b.isLeftC = "<< (b->IsParentsLeftChild() ? "T" : "F");
    o << " b.fullEdgeInfo = [";
    writeEdgeDecompInfo(o, *(b->GetFullEdgeInfoConstPtr()));
    o << "]";
    if (b->GetFullEdgeInfoConstPtr() == b->GetActiveEdgeInfoConstPtr()) {
    o << " b.activeEdgeInfo = b.fullEdgeInfo";
    }
    else {
        o << " b.activeEdgeInfo = [";
        writeEdgeDecompInfo(o, *b->GetActiveEdgeInfoConstPtr());
        o << "]";
    }
}

inline void debugNdPrint(std::ostream & o, const NxsSimpleNode & nd) {
    o << "address = " << (long) &nd;
    o << " taxon_ind = " << nd.GetTaxonIndex() << " [";
    NdBlob * blob = (NdBlob*) nd.scratch;
    if (blob) {
        debugBlobPrint(o, blob);
    }
    o << ']';
    const NxsSimpleNode *par = nd.GetParent();
    if (par) 
        o << " par = " << par->GetTaxonIndex();
    else
        o << " no parent";
    const NxsSimpleNode *lc = nd.GetFirstChild();
    if (lc) {
        o << " left = " << lc->GetTaxonIndex();
        const NxsSimpleNode *rc = lc->GetNextSib();
        if (rc) {
            o << " right = " << rc->GetTaxonIndex();
        }
        else {
            o << " no right.";
        }
    }
    else {
        o << " no left.";
    }
}

// checks that the active dir fields and num active leaves are consistently 
//  marked for the subtree that is rooted at `nd`
bool debugCheckActiveFlags(const NxsSimpleNode * nd, long numActiveAbove) {
    assert(nd);
    NdBlob * b = (NdBlob*)nd->scratch;
    if (nd->IsTip()) {
        return true;
    }
    if (numActiveAbove > 0 && b->GetNumActiveLeavesAbove() != numActiveAbove) {
        std::cerr << "debugCheckActiveFlags failure: numActiveAbove=" << numActiveAbove << " node info:\n "; debugNdPrint(std::cerr, *nd); std::cerr << '\n';
        assert(b->GetNumActiveLeavesAbove() == numActiveAbove);
    }
    int upperDir =  (int)b->GetActiveAndFocalLeafDir() & (int)LEFT_OR_RIGHT;
    if (upperDir == (int) LEFT_OR_RIGHT) {
        const NxsSimpleNode * ln = LeftChild(*nd);
        assert(ln);
        NdBlob * lb = (NdBlob*)ln->scratch;
        const NxsSimpleNode * rn = RightChild(*nd);
        assert(rn);
        NdBlob * rb = (NdBlob*)rn->scratch;
        if (b->GetNumActiveLeavesAbove() != (lb->GetNumActiveLeavesAbove() + rb->GetNumActiveLeavesAbove())) {
            std::cerr << "debugCheckActiveFlags failure: numActiveAbove=" << numActiveAbove << " node info:\n "; debugNdPrint(std::cerr, *nd); std::cerr << '\n';
            std::cerr << "                      LeftChild: numActiveAbove=" << lb->GetNumActiveLeavesAbove() << " node info:\n "; debugNdPrint(std::cerr, *ln); std::cerr << '\n';
            std::cerr << "                      RightChild: numActiveAbove=" << rb->GetNumActiveLeavesAbove() << " node info:\n "; debugNdPrint(std::cerr, *rn); std::cerr << '\n';
            assert(b->GetNumActiveLeavesAbove() == (lb->GetNumActiveLeavesAbove() + rb->GetNumActiveLeavesAbove()));
        }
        debugCheckActiveFlags(ln, lb->GetNumActiveLeavesAbove());
        debugCheckActiveFlags(rn, rb->GetNumActiveLeavesAbove());
    }
    else if (upperDir == (int) LEFT_DIR_BIT) {
        const NxsSimpleNode * ln = LeftChild(*nd);
        assert(ln);
        NdBlob * lb = (NdBlob*)ln->scratch;
        if (b->GetNumActiveLeavesAbove() != lb->GetNumActiveLeavesAbove()) {
            std::cerr << "debugCheckActiveFlags failure: numActiveAbove=" << numActiveAbove << " node info:\n "; debugNdPrint(std::cerr, *nd); std::cerr << '\n';
            std::cerr << "                      LeftChild: numActiveAbove=" << lb->GetNumActiveLeavesAbove() << " node info:\n "; debugNdPrint(std::cerr, *ln); std::cerr << '\n';
            assert(b->GetNumActiveLeavesAbove() == lb->GetNumActiveLeavesAbove());
        }
        debugCheckActiveFlags(ln, lb->GetNumActiveLeavesAbove());
    }
    else if (upperDir == (int) LEFT_DIR_BIT) {
        const NxsSimpleNode * ln = RightChild(*nd);
        assert(ln);
        NdBlob * lb = (NdBlob*)ln->scratch;
        if (b->GetNumActiveLeavesAbove() != lb->GetNumActiveLeavesAbove()) {
            std::cerr << "debugCheckActiveFlags failure: numActiveAbove=" << numActiveAbove << " node info:\n "; debugNdPrint(std::cerr, *nd); std::cerr << '\n';
            std::cerr << "                      RightChild: numActiveAbove=" << lb->GetNumActiveLeavesAbove() << " node info:\n "; debugNdPrint(std::cerr, *ln); std::cerr << '\n';
            assert(b->GetNumActiveLeavesAbove() == lb->GetNumActiveLeavesAbove());
        }
        debugCheckActiveFlags(ln, lb->GetNumActiveLeavesAbove());
    }
    else {
        std::cerr << "debugCheckActiveFlags failure. No active children. numActiveAbove=" << numActiveAbove << " node info:\n "; debugNdPrint(std::cerr, *nd); std::cerr << '\n';
        assert(false);
    }
    return true;
}


void debugPrintTree(std::ostream &o, const NxsSimpleNode *nd, std::string indentation) {
    o << indentation;
    debugNdPrint(o, *nd);
    o << '\n';
    if (!nd->IsTip()) {
        std::string n = indentation;
        n.append(4, ' ');
        const NxsSimpleNode * child = LeftChild(*nd);
        debugPrintTree(o, child, n);
        child = RightChild(*nd);
        debugPrintTree(o, child, n);
    }
}
void debugPrintTree(std::ostream &o, const NxsSimpleNode *nd) {
    std::string mt;
    debugPrintTree(o, nd, mt);
}

unsigned gMaxNumLeafPaths = 0;

// user-controlled options as global variables...
bool gVerbose = false;
long gLeafSetIntersectionSize = 48; // we demand that the leaf set intersection size be exactly THIS_NODE_DIR
long gFloorHalfIntersectionSize = 24; // we shoot for this many above and below the centroid

long gMaxSubProblemSize = 250;
bool gUseEdgeLengths = false;
bool gIntercalate = false;

std::ostream * gOutputStream = &(std::cout);



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


inline long CentroidScore(const NxsSimpleNode &nd,  long numActive) {
    NdBlob * b = (NdBlob *)nd.scratch;
    assert(b);
    const long above = b->GetNumActiveLeavesAbove();
    const long below = numActive - above;

    long sc = (above > below ? above - below : below - above);
    if (gDebugging) {
        std::cerr << "CentroidScore(nd={";
        debugNdPrint(std::cerr, nd);
        std::cerr << "}, numActive=" << numActive << ") = "<< sc << std::endl;
    }
    return sc;
}
inline void mergeShortLeavesLists(const NxsSimpleNode &nd, TreeSweepDirection dir, const NxsSimpleNode &o, const NxsSimpleNode &o2) {
    if (gIntercalate) {
        assert(false);
        throw 5;
    }
    else {
        EdgeDecompInfo * edgeInfo = ((NdBlob*)nd.scratch)->GetActiveEdgeInfoPtr();
        NdBlob* oBlob = (NdBlob*)(o.scratch);
        NdBlob* o2Blob = (NdBlob*)(o2.scratch);
        assert(edgeInfo);
        if (dir == LEFT_OR_RIGHT) {
            //std::cerr << "mergeShortLeavesLists nd = " << (long) &nd << " ndblob = " << (long) nd.scratch << " ndBlob->GetActiveEdgeInfoPtr() = " << (long)edgeInfo << " edgeInfo->closestLeavesAbove.size() = " << edgeInfo->closestLeavesAbove.size() << std::endl;
            //std::cerr << "mergeShortLeavesLists o = " << (long) &o << " oblob = " << (long) oBlob << " oBlob->GetActiveEdgeInfoPtr() = " << (long)oBlob->GetActiveEdgeInfoPtr() << " oBlob->GetActiveEdgeInfoPtr()->closestLeavesAbove.size() = " << oBlob->GetActiveEdgeInfoPtr()->closestLeavesAbove.size() << std::endl;
            //std::cerr << "mergeShortLeavesLists o22 = " << (long) &o2 << " o2blob = " << (long) o2Blob << " o2Blob->GetActiveEdgeInfoPtr() = " << (long)oBlob->GetActiveEdgeInfoPtr() << " oBlob->GetActiveEdgeInfoPtr()->closestLeavesAbove.size() = " << o2Blob->GetActiveEdgeInfoPtr()->closestLeavesAbove.size() << std::endl;
            LPECollection &peList = edgeInfo->closestLeavesAbove;
            assert(peList.empty());
            mergePathElementLists(peList, LEFT_DIR, oBlob->GetActiveEdgeInfoPtr()->closestLeavesAbove,
                                          RIGHT_DIR, o2Blob->GetActiveEdgeInfoPtr()->closestLeavesAbove);
        }
        else {
        
            EdgeDecompInfo * parInfo = o2Blob->GetActiveEdgeInfoPtr();
            EdgeDecompInfo * oInfo = oBlob->GetActiveEdgeInfoPtr();
            const LPECollection & parSource = parInfo->closestLeavesBelow;
            LPECollection &peList = edgeInfo->closestLeavesBelow;
            TreeDirection od = (dir == LEFT_OR_BELOW ? LEFT_DIR : RIGHT_DIR);
            mergePathElementLists(peList, od, oInfo->closestLeavesAbove, BELOW_DIR, parSource);
        }
    }
}

void nonRecursiveFlagActivePathDown(const NxsSimpleNode * currAnc,
                                  const LPECollection * newActiveLeaves,
                                  const NxsSimpleNode * & nextNd,
                                  const LPECollection * & nextLPEC,
                                  const NxsSimpleNode * & futureNd,
                                  const LPECollection * & futureLPEC,
                                  TreeSweepDirection * futureBitToAdd,
                                  TreeSweepDirection bitToAdd) {
    if (gDebugging) {
        std::cerr << "nonRecursiveFlagActivePathDown(currAnc=" << (long)currAnc->GetTaxonIndex() << ", newActiveLeaves=" << (long) newActiveLeaves << ')' << std::endl;
        if (currAnc->GetFirstChild()) {
            std::cerr << "LeftChild(*currAnc)=" << (long)currAnc->GetFirstChild()->GetTaxonIndex() << std::endl;
            if (RightChild(*currAnc)) {
                std::cerr << "RightChild(*currAnc)=" << (long)RightChild(*currAnc)->GetTaxonIndex()  << std::endl;
            }
        }
    }
    assert(currAnc);
    NdBlob * currBlob = (NdBlob*)(currAnc->scratch);
    assert(currBlob);

    currBlob->lpeScratch1.clear();
    currBlob->lpeScratch2.clear();
    nextNd = 0L;
    nextLPEC = 0L;
    futureNd = 0L;
    futureLPEC = 0L;
    assert(!currAnc->IsTip());
    bool hadLeft = false;
    bool hadRight = false;
    if (currAnc == gRoot) {
        return;
    }
    else {
        if (newActiveLeaves == 0L) {
            LPECollectionConstIt nalIt = currBlob->GetFullEdgeInfoConstPtr()->closestLeavesBelow.begin();
            currBlob->SetActiveEdgeInfoPtr(currBlob->GetFullEdgeInfoPtr());
            currBlob->SetNumActiveLeavesAbove(gNumLeaves - currBlob->numLeavesAboveEdge);
            for (; nalIt != currBlob->GetFullEdgeInfoConstPtr()->closestLeavesBelow.end(); ++nalIt) {
                const LeafPathElement & el = *nalIt;
                if (el.dirToNext == LEFT_DIR) {
                    currBlob->lpeScratch1.push_back(el);
                    hadLeft = true;
                    assert(!hadRight);
                }
                else if (el.dirToNext == RIGHT_DIR) {
                    currBlob->lpeScratch1.push_back(el);
                    hadRight = true;
                    assert(!hadLeft);
                }
                else {
                    assert(el.dirToNext == BELOW_DIR);
                    currBlob->lpeScratch2.push_back(el);
                }
            }
    
        }
        else {
            LPECollectionConstIt nalIt = newActiveLeaves->begin();
            currBlob->SetActiveEdgeInfoPtr(new EdgeDecompInfo());
            currBlob->SetNumActiveLeavesAbove(currBlob->numLeavesAboveEdge - newActiveLeaves->size());
            for (; nalIt != newActiveLeaves->end(); ++nalIt) {
                int indInCurr = nalIt->indexInNext;
                const LeafPathElement & el = currBlob->GetFullEdgeInfoConstPtr()->closestLeavesBelow.at(indInCurr);
                if (gDebugging) {
                    std::cerr << "nalIt->indexInNext = " << indInCurr << " el.leaf=" << el.leaf->GetTaxonIndex() << " nalIt->leaf=" << (long) nalIt->leaf->GetTaxonIndex() << std::endl;
                }
                assert(el.leaf == nalIt->leaf);
                if (el.dirToNext == LEFT_DIR) {
                    currBlob->lpeScratch1.push_back(el);
                    hadLeft = true;
                    assert(!hadRight);
                }
                else if (el.dirToNext == RIGHT_DIR) {
                    currBlob->lpeScratch1.push_back(el);
                    hadRight = true;
                    assert(!hadLeft);
                }
                else {
                    assert(el.dirToNext == BELOW_DIR);
                    currBlob->lpeScratch2.push_back(el);
                }
            }
        }
        if (currBlob->lpeScratch1.empty()) {
            if (currBlob->lpeScratch2.empty()) {
                currBlob->SetActiveLeafDir(bitToAdd);
                return;
            }
            futureNd = Parent(*currAnc);
            futureLPEC = &(currBlob->lpeScratch2);
        }
        else if (currBlob->lpeScratch2.empty()) {
            if (currBlob->IsParentsLeftChild()) {
                *futureBitToAdd = RIGHT_DIR_BIT;
                nextNd = RightChild(*Parent(*currAnc));
            }
            else {
                *futureBitToAdd = LEFT_DIR_BIT;
                nextNd = LeftChild(*Parent(*currAnc));
            }
            nextLPEC = &(currBlob->lpeScratch1);
            }
        else {
            if (currBlob->IsParentsLeftChild()) {
                *futureBitToAdd = RIGHT_DIR_BIT;
                nextNd = RightChild(*Parent(*currAnc));
            }
            else {
                *futureBitToAdd = LEFT_DIR_BIT;
                nextNd = LeftChild(*Parent(*currAnc));
            }
            nextLPEC = &(currBlob->lpeScratch1);
            futureNd = Parent(*currAnc);
            futureLPEC = &(currBlob->lpeScratch2);
        }
        currBlob->SetActiveLeafDir(TreeSweepDirection(BELOW_DIR_BIT | bitToAdd));
    }
}

void nonRecursiveFlagActivePathUp(const NxsSimpleNode * currAnc,
                                  const LPECollection * newActiveLeaves,
                                  const NxsSimpleNode * & nextNd,
                                  const LPECollection * & nextLPEC,
                                  const NxsSimpleNode * & futureNd,
                                  const LPECollection * & futureLPEC) {
    if (gDebugging) {
        std::cerr << "nonRecursiveFlagActivePathUp(currAnc=" << (long)currAnc->GetTaxonIndex() << ", newActiveLeaves=" << (long) newActiveLeaves << ')' << std::endl;
        if (currAnc->GetFirstChild()) {
            std::cerr << "LeftChild(*currAnc)=" << (long)currAnc->GetFirstChild()->GetTaxonIndex() << std::endl;
            std::cerr << "RightChild(*currAnc)=" << (long)RightChild(*currAnc)->GetTaxonIndex()  << std::endl;
        }
    }
    assert(currAnc);
    NdBlob * currBlob = (NdBlob*)(currAnc->scratch);
    assert(currBlob);

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
        LPECollectionConstIt nalIt = currBlob->GetFullEdgeInfoConstPtr()->closestLeavesAbove.begin();
        currBlob->SetActiveEdgeInfoPtr(currBlob->GetFullEdgeInfoPtr());
        currBlob->SetNumActiveLeavesAbove(currBlob->numLeavesAboveEdge);
        for (; nalIt != currBlob->GetFullEdgeInfoConstPtr()->closestLeavesAbove.end(); ++nalIt) {
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
        currBlob->SetActiveEdgeInfoPtr(new EdgeDecompInfo());
        currBlob->SetNumActiveLeavesAbove(newActiveLeaves->size());
        for (; nalIt != newActiveLeaves->end(); ++nalIt) {
            int indInCurr = nalIt->indexInNext;
            const LeafPathElement & el = currBlob->GetFullEdgeInfoConstPtr()->closestLeavesAbove.at(indInCurr);
            if (gDebugging) {
                std::cerr << "nalIt->indexInNext = " << indInCurr << " el.leaf=" << el.leaf->GetTaxonIndex() << " nalIt->leaf=" << nalIt->leaf->GetTaxonIndex() << std::endl;
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
        currBlob->SetActiveLeafDir(RIGHT_DIR_BIT);
        futureNd = RightChild(*currAnc);
        futureLPEC = &(currBlob->lpeScratch2);
    }
    else if (currBlob->lpeScratch2.empty()) {
        currBlob->SetActiveLeafDir(LEFT_DIR_BIT);
        nextNd = LeftChild(*currAnc);
        nextLPEC = &(currBlob->lpeScratch1);
        }
    else {
        currBlob->SetActiveLeafDir(LEFT_OR_RIGHT);
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
/// Sets focalEdgeDir for nodes along the path *EXCEPT* for currAnc
void flagActivePathUp(const NxsSimpleNode *currAnc, const LPECollection *newActiveLeaves) {
    if (gDebugging) {
        std::cerr << "flagActivePathUp(currAnc=" << (long) currAnc->GetTaxonIndex() << ", newActiveLeaves=" <<(long)newActiveLeaves << ")" << std::endl;
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

            NdBlob * futureBlob = (NdBlob *)futureNd->scratch;
            assert(futureBlob);
            futureBlob->SetFocalEdgeDir(BELOW_DIR);
           
        }
        if (nextNd != 0) {
            assert(nextLPEC != 0);
            ndStack.push(nextNd);
            lpecStack.push(nextLPEC);

            NdBlob * nextBlob = (NdBlob *)nextNd->scratch;
            assert(nextBlob);
            nextBlob->SetFocalEdgeDir(BELOW_DIR);
        }
    }
}


/// Sets focalEdgeDir for nodes along the path *EXCEPT* for currAnc
/// \returns the "active Root"
const NxsSimpleNode * flagActivePathDown(const NxsSimpleNode *currAnc, const LPECollection *newActiveLeaves, TreeSweepDirection dirOfOtherActive) {
    if (gDebugging) {
        std::cerr << "flagActivePathDown(currAnc=" << (long) currAnc << ", newActiveLeaves=" <<(long)newActiveLeaves << ")" << std::endl;
        if (newActiveLeaves) {
            std::cerr << " leaf paths:\n";
            for (LPECollectionConstIt i = newActiveLeaves->begin(); i != newActiveLeaves->end(); ++i) {
                std::cerr << "  " << *i << '\n';
            }
            std::cerr << std::endl;
        }
    }
    std::stack<const NxsSimpleNode *> upStack;
    std::stack<const LPECollection *> upLPECStack;
    std::stack<const NxsSimpleNode *> downStack;
    std::stack<const LPECollection *> downLPECStack;
    std::stack<TreeSweepDirection> downBitStack;
    const NxsSimpleNode * actRoot = currAnc;
    TreeSweepDirection futureBitToAdd = dirOfOtherActive;
    downStack.push(currAnc);
    downLPECStack.push(newActiveLeaves);
    downBitStack.push(futureBitToAdd);
    

    for (;;) {
        if (upStack.empty()) { 
            if (downStack.empty())
                break;
            const NxsSimpleNode * currNd = downStack.top();
            downStack.pop();
            assert(!downLPECStack.empty());
            assert(!downBitStack.empty());
            const LPECollection * currLPEC = downLPECStack.top();
            downLPECStack.pop();
            const TreeSweepDirection currBitToAdd = downBitStack.top();
            downBitStack.pop();
    
            const NxsSimpleNode * nextNd = 0L;
            const NxsSimpleNode * futureNd = 0L;
            const LPECollection * nextLPEC = 0L;
            const LPECollection * futureLPEC = 0L;
            nonRecursiveFlagActivePathDown(currNd, currLPEC, nextNd, nextLPEC, futureNd, futureLPEC, &futureBitToAdd, currBitToAdd);
            if (futureNd != 0) {
                actRoot = futureNd;
                assert(futureLPEC != 0);
                downStack.push(futureNd);
                downLPECStack.push(futureLPEC);
                downBitStack.push(futureBitToAdd);

                NdBlob * futureBlob = (NdBlob *)futureNd->scratch;
                assert(futureBlob);
                NdBlob * currBlob = (NdBlob *)currNd->scratch;
                assert(currBlob);
                futureBlob->SetFocalEdgeDir((currBlob->IsParentsLeftChild() ? LEFT_DIR : RIGHT_DIR));
            }
            if (nextNd != 0) {
                assert(nextLPEC != 0);
                upStack.push(nextNd);
                upLPECStack.push(nextLPEC);
                NdBlob * nextBlob = (NdBlob *) nextNd->scratch;
                assert(nextBlob);
                nextBlob->SetFocalEdgeDir(BELOW_DIR);
            }
        }
        else {
            const NxsSimpleNode * currNd = upStack.top();
            upStack.pop();
            assert(!upLPECStack.empty());
            const LPECollection * currLPEC = upLPECStack.top();
            upLPECStack.pop();
    
            const NxsSimpleNode * nextNd = 0L;
            const NxsSimpleNode * futureNd = 0L;
            const LPECollection * nextLPEC = 0L;
            const LPECollection * futureLPEC = 0L;
            nonRecursiveFlagActivePathUp(currNd, currLPEC, nextNd, nextLPEC, futureNd, futureLPEC);
            if (futureNd != 0) {
                assert(futureLPEC != 0);
                upStack.push(futureNd);
                upLPECStack.push(futureLPEC);
            }
            if (nextNd != 0) {
                assert(nextLPEC != 0);
                upStack.push(nextNd);
                upLPECStack.push(nextLPEC);
            }
        }
    }
    return actRoot;
}

void writeNewickOfActive(std::ostream & out, const NxsSimpleNode *activeRoot) {
    assert(activeRoot);
    typedef std::pair<const NxsSimpleNode *, std::stack<char> > NdSuffixPair;
    std::stack<NdSuffixPair> ndStack;
    std::stack<char> currSuffix;
    char currPrefix = '\0';
    for (;;) {
        NdBlob * currBlob = (NdBlob *)(activeRoot->scratch);
        //if (gDebugging)
        //    std::cerr << "activeRoot = " << (long) activeRoot << ", activeRoot->GetTaxonIndex()=" << activeRoot->GetTaxonIndex() << ", currBlob->activeLeafDir="  << currBlob->activeLeafDir << std::endl;
        if (activeRoot->IsTip()) {
            out << (activeRoot->GetTaxonIndex() + 1);
            activeRoot = 0L;
        }
        else {
            TreeSweepDirection fullActiveLeafDir = currBlob->GetActiveAndFocalLeafDir();
            //if (gDebugging)
            //    std::cerr << " fullActiveLeafDir = " << fullActiveLeafDir << std::endl;
            const bool leftActive = fullActiveLeafDir & LEFT_DIR_BIT;
            const bool rightActive = fullActiveLeafDir & RIGHT_DIR_BIT;
            if (leftActive && rightActive) {
                if (currPrefix != '\0') {
                    currSuffix.push(currPrefix);
                }
                if (gDebugging)
                    out << "[i" << activeRoot->GetTaxonIndex() << ']';
                currSuffix.push(')');
                ndStack.push(NdSuffixPair(RightChild(*activeRoot), currSuffix));
                currSuffix = std::stack<char>();
                currPrefix = ',';
                out << '(';
                activeRoot = LeftChild(*activeRoot);
            }
            else if (leftActive) {
                activeRoot = LeftChild(*activeRoot);
                assert(activeRoot);
            }
            else {
                assert(rightActive);
                activeRoot = RightChild(*activeRoot);
            }
        }        
        if (activeRoot == 0L) {
            if (currPrefix != '\0') {
                out << currPrefix;
                currPrefix = '\0';
            }
            while (!currSuffix.empty()) {                
                out << currSuffix.top();
                currSuffix.pop();
            }
            if (ndStack.empty()) {
                assert(currPrefix == '\0');
                break;
            }
            else {
                const NdSuffixPair & p = ndStack.top();
                activeRoot = p.first;
                currSuffix = p.second;
                ndStack.pop();
            }
        }
    }
}

const NxsSimpleNode * findCentroidForActiveLeaves(const NxsSimpleNode *nd, long totalNumActiveLeaves, long ndScore, long * retScore) {
    assert(!nd->IsTip());
    NdBlob * currBlob = (NdBlob *)(nd->scratch);
    *retScore = ndScore;
    const NxsSimpleNode * possibleReturn;
    TreeSweepDirection fullActiveLeafDir = currBlob->GetActiveAndFocalLeafDir();
    if (!nd->IsTip()) {
        const NxsSimpleNode * child;
        long childScore;
        if (fullActiveLeafDir & LEFT_DIR_BIT) {
            if (gDebugging)
                std::cerr << "Calculating LeftChild centroid score... ";
            child = LeftChild(*nd);
            childScore = CentroidScore(*child, totalNumActiveLeaves);
            if (childScore <= ndScore) {
                possibleReturn = findCentroidForActiveLeaves(child, totalNumActiveLeaves, childScore, retScore);
                if (*retScore < ndScore) {
                    if (gDebugging)
                        std::cerr << "Moving Left from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << *retScore << '\n';
                    return possibleReturn;
                }
            }
        }
        if (fullActiveLeafDir & RIGHT_DIR_BIT) {
            if (gDebugging)
                std::cerr << "Calculating RightChild centroid score... ";
            child = RightChild(*nd);
            childScore = CentroidScore(*child, totalNumActiveLeaves);
            if (childScore <= ndScore) {
                possibleReturn = findCentroidForActiveLeaves(child, totalNumActiveLeaves, childScore, retScore);
                if (*retScore < ndScore) {
                    if (gDebugging)
                        std::cerr << "Moving Right from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << *retScore << '\n';
                    return possibleReturn;
                }
            }
        }
    }
    if (fullActiveLeafDir & BELOW_DIR_BIT ) {
        const NxsSimpleNode * par = Parent(*nd);
        assert(par);
        if (gDebugging)
            std::cerr << "Calculating Parent centroid score... ";
        long parScore = CentroidScore(*par, totalNumActiveLeaves);
        if (parScore <= ndScore) {
            possibleReturn =  findCentroidForActiveLeaves(nd, totalNumActiveLeaves, parScore, retScore);
            if (*retScore < ndScore) {
                if (gDebugging)
                    std::cerr << "Moving Down from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << *retScore << '\n';
                return possibleReturn;
            }
        }
    }
    if (gDebugging)
        std::cerr << "Centroid Edge found nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ")\n";
    return nd;
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
void decomposeAroundCentroidChild(std::vector<const NxsSimpleNode *> &preorderTraversal,
                                  const NxsSimpleNode *topCentroidChild,
                                  const NxsSimpleNode *bottomCentroidChild,
                                  long numActiveLeaves,
                                  NxsString namePrefix) {
    ////////////////////////////////////////////////////////////////////////////
    // Step one is to get identify leftSubtreeRoot, rightSubtreeRoot, sibSubtreeRoot
    //  and gpSubtreeRoot (see diagrams above).
    ////////////////////////////////////////////////////////////////////////////
    assert(topCentroidChild);
    assert(!topCentroidChild->IsTip());
    assert(bottomCentroidChild);
    assert(!bottomCentroidChild->IsTip());
    const NxsSimpleNode * leftSubtreeRoot = LeftChild(*topCentroidChild);
    assert(leftSubtreeRoot);
    const NxsSimpleNode * rightSubtreeRoot = RightChild(*topCentroidChild);
    assert(rightSubtreeRoot);
    const NxsSimpleNode * parSubtreeRoot = Parent(*bottomCentroidChild);
    assert(parSubtreeRoot);
    assert(!parSubtreeRoot->IsTip());
    const NxsSimpleNode * parSubtreeRootlc = LeftChild(*parSubtreeRoot);
    const NxsSimpleNode * sibSubtreeRoot = 0L;
    if (parSubtreeRootlc == bottomCentroidChild) {
        sibSubtreeRoot = parSubtreeRootlc->GetNextSib();
        assert(sibSubtreeRoot);
        assert(sibSubtreeRoot->GetNextSib() == 0L);
    }
    else {
        assert(parSubtreeRootlc->GetNextSib() == bottomCentroidChild);
        assert(bottomCentroidChild->GetNextSib() == 0L);
        sibSubtreeRoot = parSubtreeRootlc;
    }
    assert(sibSubtreeRoot);

/*    const NxsSimpleNode * gpSubtreeRoot = Parent(*parSubtreeRoot);
    assert(gpSubtreeRoot);
    if (IsRootChild(*parSubtreeRoot)) {
        gpSubtreeRoot = LeftChild(*Parent(*parSubtreeRoot));
        assert(gpSubtreeRoot->GetNextSib() == parSubtreeRoot);
        assert(parSubtreeRoot->GetNextSib() == 0L);
    }
    else {
        assert(Parent(*gpSubtreeRoot));
    }
*/
    NdBlob * centroidBlob = (NdBlob *)(topCentroidChild->scratch);
    NdBlob * parBlob = (NdBlob *)(parSubtreeRoot->scratch); 
    NdBlob * sibBlob = (NdBlob *)(sibSubtreeRoot->scratch); 
    NdBlob * leftBlob = (NdBlob *)(leftSubtreeRoot->scratch);
    NdBlob * rightBlob = (NdBlob *)(rightSubtreeRoot->scratch);
    leftBlob->SetFocalEdgeDir(BELOW_DIR);
    rightBlob->SetFocalEdgeDir(BELOW_DIR);
    sibBlob->SetFocalEdgeDir(BELOW_DIR);
    parBlob->SetFocalEdgeDir(centroidBlob->IsParentsLeftChild() ? LEFT_DIR : RIGHT_DIR);
    
    ////////////////////////////////////////////////////////////////////////////
    // Step two is to identify the common leaf set for all four decompositions.
    ////////////////////////////////////////////////////////////////////////////
    long numActLeavesAboveCentroid = centroidBlob->GetNumActiveLeavesAbove();
    long numActLeavesBelowCentroid = numActiveLeaves - numActLeavesAboveCentroid;
    long numAboveChosen, numBelowChosen, numLeftChosen, numRightChosen, numSibChosen, numParChosen;
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
    const long origNumActiveLeavesAboveLeft = leftBlob->GetNumActiveLeavesAbove();
    const long origNumActiveLeavesAboveRight = rightBlob->GetNumActiveLeavesAbove();
    const long origNumActiveLeavesAboveSib = sibBlob->GetNumActiveLeavesAbove();
    if (origNumActiveLeavesAboveLeft < floorHalfAbove) {
        numLeftChosen = origNumActiveLeavesAboveLeft;
        assert(rightBlob->GetNumActiveLeavesAbove() >= floorHalfAbove - numLeftChosen);
        numRightChosen = numAboveChosen - numLeftChosen;
    }
    else {
        if (rightBlob->GetNumActiveLeavesAbove() < floorHalfAbove) {
            numRightChosen = rightBlob->GetNumActiveLeavesAbove();
            assert(origNumActiveLeavesAboveLeft >= numAboveChosen - numRightChosen);

        }
        else {
            numRightChosen = floorHalfAbove; //\TODO should probably sort the PLE's and use their order to decide which direction gets the rounding error...
        }
        numLeftChosen = numAboveChosen - numRightChosen;
    }
    if (gDebugging)
        std::cerr << "numAboveChosen=" <<numAboveChosen << ", numBelowChosen=" << numBelowChosen << " sibBlob->GetNumActiveLeavesAbove()=" << sibBlob->GetNumActiveLeavesAbove() <<  ", numActiveLeaves=" << numActiveLeaves << ", parBlob->GetNumActiveLeavesAbove()=" << parBlob->GetNumActiveLeavesAbove() << '\n';


    long floorHalfBelow = numBelowChosen/2;
    if (sibBlob->GetNumActiveLeavesAbove() < floorHalfBelow) {
        numSibChosen = sibBlob->GetNumActiveLeavesAbove();
        assert((numActiveLeaves  - parBlob->GetNumActiveLeavesAbove()) >= numBelowChosen - numSibChosen);
        numParChosen = numBelowChosen - numSibChosen;
    }
    else {
        if ((numActiveLeaves  - parBlob->GetNumActiveLeavesAbove()) < floorHalfBelow) {
            numParChosen = (numActiveLeaves  - parBlob->GetNumActiveLeavesAbove());
            assert(sibBlob->GetNumActiveLeavesAbove() == numBelowChosen - numParChosen);
        }
        else {
            numParChosen = floorHalfBelow; //\TODO should probably sort the PLE's and use their order to decide which direction gets the rounding error...
        }
        numSibChosen = numBelowChosen - numParChosen;
    }

    assert(numAboveChosen + numBelowChosen == gLeafSetIntersectionSize);
    assert(numLeftChosen + numRightChosen == numAboveChosen);
    assert(numSibChosen + numParChosen == numBelowChosen);
    assert(numLeftChosen > 0);
    assert(numRightChosen > 0);
    assert(numSibChosen > 0);
    assert(numParChosen > 0);


    const EdgeDecompInfo * leftEDI = leftBlob->GetActiveEdgeInfoConstPtr();
    assert(leftEDI);
    const LPECollection & leftFullLPEC = leftEDI->closestLeavesAbove;
    assert(leftFullLPEC.size() >= (unsigned long)numLeftChosen);
    LPECollection::const_iterator lfLPECIt = leftFullLPEC.begin();
    LPECollection leftCommonLeafSet;
    leftCommonLeafSet.reserve(numLeftChosen);
    for (unsigned i=0; i < numLeftChosen; ++i, ++lfLPECIt) {
        leftCommonLeafSet.push_back(LeafPathElement(*lfLPECIt, LEFT_DIR, i));
    }
    assert(leftCommonLeafSet.size() == (unsigned long) numLeftChosen);

    const EdgeDecompInfo * rightEDI = rightBlob->GetActiveEdgeInfoConstPtr();
    assert(rightEDI);
    const LPECollection & rightFullLPEC = rightEDI->closestLeavesAbove;
    assert(rightFullLPEC.size() >= (unsigned long)numRightChosen);
    LPECollection::const_iterator rfLPECIt = rightFullLPEC.begin();
    LPECollection rightCommonLeafSet;
    rightCommonLeafSet.reserve(numRightChosen);
    for (unsigned i=0; i < numRightChosen; ++i, ++rfLPECIt) {
        rightCommonLeafSet.push_back(LeafPathElement(*rfLPECIt, RIGHT_DIR, i));
    }
    assert(rightCommonLeafSet.size() == (unsigned long) numRightChosen);

    const EdgeDecompInfo * sibEDI = sibBlob->GetActiveEdgeInfoConstPtr();
    assert(sibEDI);
    const LPECollection & sibFullLPEC = sibEDI->closestLeavesAbove;
    assert(sibFullLPEC.size() >= (unsigned long)numSibChosen);
    LPECollection::const_iterator sibfLPECIt = sibFullLPEC.begin();
    LPECollection sibCommonLeafSet;
    sibCommonLeafSet.reserve(numSibChosen);
    for (unsigned i=0; i < numSibChosen; ++i, ++sibfLPECIt) {
        sibCommonLeafSet.push_back(LeafPathElement(*sibfLPECIt, (sibBlob->IsParentsLeftChild() ? LEFT_DIR : RIGHT_DIR), i));
    }
    assert(sibCommonLeafSet.size() == (unsigned long)numSibChosen);

    const EdgeDecompInfo * parEDI = parBlob->GetActiveEdgeInfoConstPtr();
    assert(parEDI);
    const LPECollection & parFullLPEC = parEDI->closestLeavesBelow;
    assert(parFullLPEC.size() >= (unsigned long)numParChosen);
    LPECollection::const_iterator parfLPECIt = parFullLPEC.begin();
    LPECollection parCommonLeafSet;
    parCommonLeafSet.reserve(numParChosen);
    for (unsigned i=0; i < numParChosen; ++i, ++parfLPECIt) {
        parCommonLeafSet.push_back(LeafPathElement(*parfLPECIt, BELOW_DIR, i));
    }
    assert(parCommonLeafSet.size() == (unsigned long)numParChosen);

    ////////////////////////////////////////////////////////////////////////////
    // Step 3 - 6 activate each subproblem in turn and recurse...
    //
    ////////////////////////////////////////////////////////////////////////////
    if (gDebugging) {
        std::cerr << "Before flagging leftBlob->GetNumActiveLeavesAbove() = " << origNumActiveLeavesAboveLeft << ", currBlob->GetFullEdgeInfoConstPtr()->closestLeavesAbove.size() = " << leftBlob->GetFullEdgeInfoConstPtr()->closestLeavesAbove.size() << std::endl;
    }
    debugCheckActiveFlags(gRoot, gNumLeaves);
    flagActivePathUp(leftSubtreeRoot, 0L);
    debugCheckActiveFlags(leftSubtreeRoot, origNumActiveLeavesAboveLeft);
    flagActivePathUp(rightSubtreeRoot, &rightCommonLeafSet);
    debugCheckActiveFlags(rightSubtreeRoot, numRightChosen);
    flagActivePathUp(sibSubtreeRoot, &sibCommonLeafSet);
    debugCheckActiveFlags(sibSubtreeRoot, numSibChosen);
    TreeSweepDirection parOtherActiveDir = sibBlob->IsParentsLeftChild() ? LEFT_DIR_BIT : RIGHT_DIR_BIT;
    const NxsSimpleNode * activeRoot = flagActivePathDown(parSubtreeRoot, &parCommonLeafSet, parOtherActiveDir);
    const NxsSimpleNode * tmpNd = topCentroidChild;
    while (tmpNd != parSubtreeRoot) {
        NdBlob * tmpBlob = (NdBlob *)tmpNd->scratch;
        tmpBlob->SetNumActiveLeavesAbove(numAboveChosen);
        tmpNd = Parent(*tmpNd);
    }
    long sizeOfLeftSubproblem = numBelowChosen + numRightChosen + origNumActiveLeavesAboveLeft;

    if (gDebugging) {
        std::cerr << "numLeftChosen = " << numLeftChosen << ", numRightChosen = " << numRightChosen  << ", numSibChosen = " << numSibChosen << ",  numParChosen = " << numParChosen << std::endl;
        std::cerr << "leftBlob->GetNumActiveLeavesAbove() = " << leftBlob->GetNumActiveLeavesAbove() << ", sizeOfLeftSubproblem = " << sizeOfLeftSubproblem << std::endl;
        std::cerr << "Common leafset: left: ";
        writeLeafSet(std::cerr, leftCommonLeafSet);
        std::cerr << "\nCommon leafset: right: ";
        writeLeafSet(std::cerr, rightCommonLeafSet);
        std::cerr << "\nCommon leafset: sib: ";
        writeLeafSet(std::cerr, sibCommonLeafSet);
        std::cerr << "\nCommon leafset: parSubtreeRoot: ";
        writeLeafSet(std::cerr, parCommonLeafSet);
        std::cerr << "\n";

        std::cerr << "topCentroidChild Node: ";
        debugNdPrint(std::cerr, *topCentroidChild);
        std::cerr << "\n";

        std::cerr << "bottomCentroidChild Node: ";
        debugNdPrint(std::cerr, *bottomCentroidChild);
        std::cerr << "\n";

        std::cerr << "left Node: ";
        debugNdPrint(std::cerr, *leftSubtreeRoot);
        std::cerr << "\n";


        std::cerr << "right Node: ";
        debugNdPrint(std::cerr, *rightSubtreeRoot);
        std::cerr << "\n";


        std::cerr << "sib Node: ";
        debugNdPrint(std::cerr, *sibSubtreeRoot);
        std::cerr << "\n";


        std::cerr << "par Node: ";
        debugNdPrint(std::cerr, *parSubtreeRoot);
        std::cerr << "\n";
        debugPrintTree(std::cerr, gRoot);
        debugCheckActiveFlags(activeRoot, -1); //sizeOfLeftSubproblem);
    }
    
        
    
    NxsString n = namePrefix;
    n += ".0";
    if (gOutputStream != 0L && (sizeOfLeftSubproblem <= gMaxSubProblemSize || gEmitIntermediates)) {
        *gOutputStream << "Tree " << n << " = [&U] ";
        writeNewickOfActive(*gOutputStream, activeRoot);
        *gOutputStream << ";\n";
    }
    
    
    if (sizeOfLeftSubproblem > gMaxSubProblemSize) {
        long sc = CentroidScore(*rootLeft, sizeOfLeftSubproblem);
        long minCS;
        const NxsSimpleNode * nextCentroid = findCentroidForActiveLeaves(topCentroidChild, sizeOfLeftSubproblem, sc, &minCS);
        const NxsSimpleNode * nextBottomCentroid = nextCentroid;
        for (;;) {
            const NxsSimpleNode * parnb  = Parent(*nextBottomCentroid);
            NdBlob * currBlob = (NdBlob *)parnb->scratch;
            TreeSweepDirection fullActiveLeafDir = currBlob->GetActiveAndFocalLeafDir();
            if ((fullActiveLeafDir & LEFT_DIR_BIT) && (fullActiveLeafDir & RIGHT_DIR_BIT)) {
                break;
            }
            nextBottomCentroid = parnb;
        }

        const NxsSimpleNode * nextTopCentroid = nextCentroid;
        for (;;) {
            NdBlob * currBlob = (NdBlob *)nextTopCentroid->scratch;
            TreeSweepDirection fullActiveLeafDir = currBlob->GetActiveAndFocalLeafDir();
            if (fullActiveLeafDir  == LEFT_OR_RIGHT) {
                break;
            }
            if (fullActiveLeafDir == LEFT_DIR_BIT) {
                nextTopCentroid = LeftChild(*nextTopCentroid);
            }
            else {
                assert (fullActiveLeafDir == RIGHT_DIR_BIT);
                nextTopCentroid = RightChild(*nextTopCentroid);
            }
        }

        decomposeAroundCentroidChild(preorderTraversal, nextTopCentroid, nextBottomCentroid, sizeOfLeftSubproblem, n);
    }
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
    rootLeft = const_cast<NxsSimpleNode *>(LeftChild(*iniRoot));
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
    unsigned currInternalIndex = gNumLeaves;
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
            rlb->SetIsParentsLeftChild(true);
            rootLeft->scratch = (void *) rlb;
            rlb->numLeavesAboveEdge = 1;
            rlb->SetNumActiveLeavesAbove(1);
            EdgeDecompInfo * rlEdgeInfo = rlb->GetActiveEdgeInfoPtr();
            rlEdgeInfo->closestLeavesAbove.push_back(LeafPathElement(rootLeft, THIS_NODE_DIR));
            rlb->SetActiveLeafDir(ALL_DIR_BITS);

            std::vector<const NxsSimpleNode *>::const_reverse_iterator postIt = preorderTraversal.rbegin();

            for (; postIt != preorderTraversal.rend(); ++postIt) {
                const NxsSimpleNode * nd = *postIt;
                NdBlob * nb = new NdBlob(LeftChild(*Parent(*nd)) == nd);
                gAllocedBlobs.push_back(nb);
                nd->scratch = (void *) nb;
                nb->SetActiveLeafDir(ALL_DIR_BITS);

                if (nd->IsTip()) {
                    nb->numLeavesAboveEdge = 1;
                    nb->SetNumActiveLeavesAbove(1);
                    EdgeDecompInfo * edgeInfo = nb->GetActiveEdgeInfoPtr();
                    edgeInfo->closestLeavesAbove.push_back(LeafPathElement(nd, THIS_NODE_DIR));
                }
                else {
                    nb->numLeavesAboveEdge = LeftBlob(*nd)->numLeavesAboveEdge + RightBlob(*nd)->numLeavesAboveEdge;
                    nb->SetNumActiveLeavesAbove(nb->numLeavesAboveEdge);
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
            rb->SetIsParentsLeftChild(false);
            gRoot->scratch = (void *) rb;
            rb->numLeavesAboveEdge = LeftBlob(*gRoot)->numLeavesAboveEdge + RightBlob(*gRoot)->numLeavesAboveEdge;
            rb->SetNumActiveLeavesAbove(rb->numLeavesAboveEdge);
            rb->SetActiveLeafDir(ALL_DIR_BITS);
            gRoot->SetTaxonIndex(currInternalIndex++);
            
            mergeShortLeavesLists(*gRoot, LEFT_OR_RIGHT, *LeftChild(*gRoot), *RightChild(*gRoot));

            if (((NdBlob *)gRoot->scratch)->numLeavesAboveEdge != gNumLeaves) {
                gErrorMessage << "Expected " << gNumLeaves << " leaves at the root, but only found " << ((NdBlob *)gRoot->scratch)->numLeavesAboveEdge;
                throw NxsException(gErrorMessage);
            }

            std::vector<const NxsSimpleNode *>::const_iterator preIt = preorderTraversal.begin();
            for (; preIt != preorderTraversal.end(); ++preIt) {
                const NxsSimpleNode * nd = *preIt;
                assert(nd != 0L);
                if (!nd->IsTip()) {
                    (const_cast<NxsSimpleNode *>(nd))->SetTaxonIndex(currInternalIndex++);

                    NdBlob* ndBlob = (NdBlob*)(nd->scratch);
                    EdgeDecompInfo * edgeInfo = ndBlob->GetActiveEdgeInfoPtr();
                    assert(edgeInfo);
                    if (IsRootChild(*nd)) {
                        const NxsSimpleNode * leafChild = gRoot->GetFirstChild();
                        assert(leafChild != nd);
                        assert(leafChild->GetNextSib() == nd);
                        assert(leafChild->IsTip());
                        EdgeDecompInfo * lcEdgeInfo = ((NdBlob*)leafChild->scratch)->GetActiveEdgeInfoPtr();
                        assert(lcEdgeInfo);
                        assert(lcEdgeInfo->closestLeavesAbove.size() == 1);
                        assert(edgeInfo->closestLeavesBelow.empty());
                        edgeInfo->closestLeavesBelow.push_back(LeafPathElement(lcEdgeInfo->closestLeavesAbove[0], BELOW_DIR, 0));
                    }
                    else {
                        const NxsSimpleNode * par = Parent(*nd);
                        if (ndBlob->IsParentsLeftChild())
                            mergeShortLeavesLists(*nd, RIGHT_OR_BELOW, *RightChild(*par), *par);
                        else
                            mergeShortLeavesLists(*nd, LEFT_OR_BELOW, *LeftChild(*par), *par);
                    }
                }
                if (gDebugging) {
                    NdBlob * nb = (NdBlob *) nd->scratch;
                    std::cerr << "nd = " << (long) nd << " blob=" << (long) nb << " blob->numLeavesAboveEdge = " << nb->numLeavesAboveEdge <<  " blob->GetActiveEdgeInfoConstPtr() = " << (long) nb->GetActiveEdgeInfoConstPtr() <<  " blob->GetActiveEdgeInfoConstPtr()->closestLeavesAbove.size() = " << nb->GetActiveEdgeInfoConstPtr()->closestLeavesAbove.size() << '\n';
                    if (nd->IsTip()) {
                        std::cerr << "nd is tip, so below info not calculated" << std::endl;
                    }
                    else {
                        std::cerr << "nd = " << (long) nd << " blob=" << (long) nb << " blob->numLeavesBelowEdge = " << (gNumLeaves - nb->numLeavesAboveEdge) <<  " blob->GetActiveEdgeInfoConstPtr() = " << (long) nb->GetActiveEdgeInfoConstPtr() <<  " blob->GetActiveEdgeInfoConstPtr()->closestLeavesBelow.size() = " << nb->GetActiveEdgeInfoConstPtr()->closestLeavesBelow.size() << '\n';
                    }
                }
            }
        }

        assert(centroidChild != 0);
        NxsString namePrefix;
        namePrefix += gCurrTreeIndex;
        if (gEmitIntermediates && gOutputStream != 0L) {
            *gOutputStream << "Tree " << namePrefix << " = [&U] ";
            writeNewickOfActive(*gOutputStream, Parent(*rootLeft));
            *gOutputStream << ";\n";
        }

        if (gDebugging) {
            debugPrintTree(std::cerr, gRoot);
        }

        decomposeAroundCentroidChild(preorderTraversal, centroidChild, centroidChild, gNumLeaves, namePrefix);
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
			gEmitIntermediates = true;
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


