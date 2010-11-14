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

unsigned gMaxNumLeafPaths = 0;

// user-controlled options as global variables...
bool gVerbose = false;
long gLeafSetIntersectionSize = 48; // we demand that the leaf set intersection size be exactly THIS_NODE_DIR
long gFloorHalfIntersectionSize = 24; // we shoot for this many above and below the centroid

long gMaxSubProblemSize = 250;
bool gUseEdgeLengths = false;
bool gIntercalate = false;

std::ostream * gOutputStream = &(std::cout);


const LPECollection & GetClosestLeavesAboveNd(const NxsSimpleNode *nd);

inline const LPECollection & GetClosestLeavesAboveNd(const NxsSimpleNode *nd) {
    NdBlob * b = CurrentBlob(*nd);
    EdgeDecompInfo * edi = b->GetActiveEdgeInfoPtr();
    if (!edi->IsAboveInitialized()) {
        TreeSweepDirection d = b->GetActiveAndFocalLeafDir();
        TreeSweepDirection u = TreeSweepDirection(d & LEFT_OR_RIGHT);
        const LPECollection emptyCollection;
        if (u == LEFT_OR_RIGHT) {
            mergePathElementLists(edi->GetClosestLeavesAboveRef(), 
                                  LEFT_DIR, GetClosestLeavesAboveNd(LeftChild(*nd)),
                                  RIGHT_DIR, GetClosestLeavesAboveNd(RightChild(*nd)));
        }
        else if (u == LEFT_DIR_BIT) {
            mergePathElementLists(edi->GetClosestLeavesAboveRef(), 
                                  LEFT_DIR, GetClosestLeavesAboveNd(LeftChild(*nd)),
                                  RIGHT_DIR, emptyCollection);
        }
        else {
            assert(u == RIGHT_DIR_BIT);
            mergePathElementLists(edi->GetClosestLeavesAboveRef(), 
                                  LEFT_DIR, emptyCollection,
                                  RIGHT_DIR, GetClosestLeavesAboveNd(RightChild(*nd)));
        }
        edi->SetAboveInitialized(true);
    }

    if (gDebugging) {
        if (edi->GetCloseLeavesAbove().size() < (unsigned) b->GetNumActiveLeavesAbove() && b->GetNumActiveLeavesAbove() < gFloorHalfIntersectionSize) {
            std::cerr << "edi->GetCloseLeavesAbove().size() = " << edi->GetCloseLeavesAbove().size() << " b->GetNumActiveLeavesAbove() = " << b->GetNumActiveLeavesAbove() << std::endl;
            assert(false);
        }
    }
    return edi->GetCloseLeavesAbove();
    
}

inline const LPECollection & GetClosestLeavesBelowNd(const NxsSimpleNode *nd) {
    NdBlob * b = CurrentBlob(*nd);
    EdgeDecompInfo * edi = b->GetActiveEdgeInfoPtr();
    if (!edi->IsBelowInitialized()) {
        const NxsSimpleNode * par = Parent(*nd);
        assert(par);
        TreeSweepDirection d = CurrentBlob(*par)->GetActiveAndFocalLeafDir();
        const LPECollection emptyCollection;
        if (b->IsParentsLeftChild()) {
            TreeSweepDirection u = TreeSweepDirection(d & RIGHT_OR_BELOW);

            if (u == RIGHT_OR_BELOW) {
                mergePathElementLists(edi->GetClosestLeavesBelowRef(), 
                                      RIGHT_DIR, GetClosestLeavesAboveNd(RightChild(*par)),
                                      BELOW_DIR, GetClosestLeavesBelowNd(par));
            }
            else if (u == BELOW_DIR_BIT) {
                mergePathElementLists(edi->GetClosestLeavesBelowRef(), 
                                      RIGHT_DIR, emptyCollection,
                                      BELOW_DIR, GetClosestLeavesBelowNd(par));
            }
            else {
                assert(u == LEFT_DIR_BIT);
                mergePathElementLists(edi->GetClosestLeavesBelowRef(), 
                                      RIGHT_DIR, GetClosestLeavesAboveNd(RightChild(*par)),
                                      BELOW_DIR, emptyCollection);
            }


        }
        else {
            TreeSweepDirection u = TreeSweepDirection(d & LEFT_OR_BELOW);

            if (u == LEFT_OR_BELOW) {
                mergePathElementLists(edi->GetClosestLeavesBelowRef(), 
                                      LEFT_DIR, GetClosestLeavesAboveNd(LeftChild(*par)),
                                      BELOW_DIR, GetClosestLeavesBelowNd(par));
            }
            else if (u == BELOW_DIR_BIT) {
                mergePathElementLists(edi->GetClosestLeavesBelowRef(), 
                                      LEFT_DIR, emptyCollection,
                                      BELOW_DIR, GetClosestLeavesBelowNd(par));
            }
            else {
                assert(u == LEFT_DIR_BIT);
                mergePathElementLists(edi->GetClosestLeavesBelowRef(), 
                                      LEFT_DIR, GetClosestLeavesAboveNd(LeftChild(*par)),
                                      BELOW_DIR, emptyCollection);
            }

        }
        edi->SetBelowInitialized(true);
    }
    return edi->GetCloseLeavesBelow();
    
}



inline void NdBlob::Reset() {
    this->numLeavesAboveEdge = -1;
    assert(edgeInfoStack.size() < 2);
    if (edgeInfoStack.size() == 1) {
        delete edgeInfoStack.top();
        edgeInfoStack.pop();
    }
    assert(edgeInfoStack.empty());
    activeEdgeInfo = &(this->fullEdgeInfo);

    this->edgeInfoStack.push(activeEdgeInfo);
    this->focalEdgeDirStack.push(this->focalEdgeDir);
    this->activeLeafDirStack.push(this->activeLeafDir);
}

inline void NdBlob::popBlob() {
    if (gDebugging) {
        std::cerr << "popped before ";
        debugBlobPrint(std::cerr, this);
        std::cerr << '\n';
    }
    if (this->activeEdgeInfo != &(this->fullEdgeInfo)) {
        delete this->activeEdgeInfo;
    }
    assert(!this->edgeInfoStack.empty());
    this->activeEdgeInfo = this->edgeInfoStack.top();
    this->edgeInfoStack.pop();

    assert(!this->activeLeafDirStack.empty());
    this->activeLeafDir  = this->activeLeafDirStack.top();
    this->activeLeafDirStack.pop();
    
    assert(!this->focalEdgeDirStack.empty());
    this->focalEdgeDir = this->focalEdgeDirStack.top();
    this->focalEdgeDirStack.pop();
    
    assert(!this->numActiveAboveStack.empty());
    this->numActiveLeavesAboveEdge = this->numActiveAboveStack.top();
    this->numActiveAboveStack.pop();
    if (gDebugging) {
        std::cerr << "popped after " ;
        debugBlobPrint(std::cerr, this);
        std::cerr << '\n';

    }
}



inline void debugNdPrint(std::ostream & o, const NxsSimpleNode & nd) {
    o << "address = " << (long) &nd;
    o << " taxon_ind = " << nd.GetTaxonIndex() << " [";
    NdBlob * blob = CurrentBlob(nd);
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
    NdBlob * b = CurrentBlob(*nd);
    if (nd->IsTip()) {
        return true;
    }
    const NxsSimpleNode * ln = LeftChild(*nd);
    const NxsSimpleNode * rn = RightChild(*nd);
    assert(ln);
    assert(rn);
    NdBlob * lb = CurrentBlob(*ln);
    if (nd == gRoot) {
        if (lb->GetNumActiveLeavesAbove() == 1)
            return debugCheckActiveFlags(ln, 1) && debugCheckActiveFlags(rn, numActiveAbove - 1);
        return debugCheckActiveFlags(ln, 0) && debugCheckActiveFlags(rn, numActiveAbove);
    }
    NdBlob * rb = CurrentBlob(*rn);
    const NxsSimpleNode * sn = SibNode(*nd);
    NdBlob * sb = CurrentBlob(*sn);
    const NxsSimpleNode * pn = Parent(*nd);
    NdBlob * pb = CurrentBlob(*pn);
    
    const EdgeDecompInfo * edi = b->GetFullEdgeInfoConstPtr();
    const LPECollection & lpec = edi->GetCloseLeavesAbove();
    for (LPECollectionConstIt i = lpec.begin(); i != lpec.end(); ++i) {
        const LeafPathElement & lpe = *i;
        if (lpe.dirToNext == LEFT_DIR) {
            const LPECollection & ld = lb->GetFullEdgeInfoConstPtr()->GetCloseLeavesAbove();
            assert(lpe.GetLeaf() == ld[lpe.indexInNext].GetLeaf());
        }
        else {
            assert(lpe.dirToNext == RIGHT_DIR);
            assert(lpe.GetLeaf() == rb->GetFullEdgeInfoConstPtr()->GetCloseLeavesAbove().at(lpe.indexInNext).GetLeaf());
        }
    }
    bool hadLeft = false;
    bool hadRight = false;
    for (LPECollectionConstIt i = edi->GetCloseLeavesBelow().begin(); i != edi->GetCloseLeavesBelow().end(); ++i) {
        const LeafPathElement & lpe = *i;
        if (lpe.dirToNext == LEFT_DIR) {
            assert(lpe.GetLeaf() == sb->GetFullEdgeInfoConstPtr()->GetCloseLeavesAbove().at(lpe.indexInNext).GetLeaf());
            hadLeft = true;
            assert(!hadRight);
        }
        else if (lpe.dirToNext == RIGHT_DIR) {
            assert(lpe.GetLeaf() == sb->GetFullEdgeInfoConstPtr()->GetCloseLeavesAbove().at(lpe.indexInNext).GetLeaf());
            hadRight = true;
            assert(!hadLeft);
        }
        else {
            assert(lpe.GetLeaf() == pb->GetFullEdgeInfoConstPtr()->GetCloseLeavesBelow().at(lpe.indexInNext).GetLeaf());
        }
    }


    if (numActiveAbove > 0 && b->GetNumActiveLeavesAbove() != numActiveAbove) {
        std::cerr << "debugCheckActiveFlags failure: numActiveAbove=" << numActiveAbove << " node info:\n "; debugNdPrint(std::cerr, *nd); std::cerr << '\n';
        assert(b->GetNumActiveLeavesAbove() == numActiveAbove);
    }
    int upperDir =  (int)b->GetActiveAndFocalLeafDir() & (int)LEFT_OR_RIGHT;
    if (upperDir == (int) LEFT_OR_RIGHT) {
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
        NdBlob * lb = CurrentBlob(*ln);
        if (b->GetNumActiveLeavesAbove() != lb->GetNumActiveLeavesAbove()) {
            std::cerr << "debugCheckActiveFlags failure: numActiveAbove=" << numActiveAbove << " node info:\n "; debugNdPrint(std::cerr, *nd); std::cerr << '\n';
            std::cerr << "                      LeftChild: numActiveAbove=" << lb->GetNumActiveLeavesAbove() << " node info:\n "; debugNdPrint(std::cerr, *ln); std::cerr << '\n';
            assert(b->GetNumActiveLeavesAbove() == lb->GetNumActiveLeavesAbove());
        }
        debugCheckActiveFlags(ln, lb->GetNumActiveLeavesAbove());
    }
    else if (upperDir == (int) RIGHT_DIR_BIT) {
        if (b->GetNumActiveLeavesAbove() != rb->GetNumActiveLeavesAbove()) {
            std::cerr << "debugCheckActiveFlags failure: numActiveAbove=" << numActiveAbove << " node info:\n "; debugNdPrint(std::cerr, *nd); std::cerr << '\n';
            std::cerr << "                      RightChild: numActiveAbove=" << rb->GetNumActiveLeavesAbove() << " node info:\n "; debugNdPrint(std::cerr, *rn); std::cerr << '\n';
            assert(b->GetNumActiveLeavesAbove() == rb->GetNumActiveLeavesAbove());
        }
        debugCheckActiveFlags(ln, lb->GetNumActiveLeavesAbove());
    }
    else {
        std::cerr << "debugCheckActiveFlags failure. No active children. upperDir=" << TreeSweepDirection(upperDir) << " node info:\n "; debugNdPrint(std::cerr, *nd); std::cerr << '\n';
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




void mergePathElementLists(LPECollection &peList,
                           TreeDirection firDir,
                           const LPECollection & firSource,
                           TreeDirection secDir,
                           const LPECollection & secSource
                           ) {
    //std::cerr << "mergePathElementLists insizes = " << peList.size() << ", " << firSource.size() << ", " << secSource.size()   << std::endl;
    peList.reserve(gMaxNumLeafPaths);
    assert(!gIntercalate);

    if (gDebugging) {
        std::cerr << "mergePathElementLists(peList, firDir=" << firDir << ", firSource=";
        writeLeafPathElementVector(std::cerr, firSource);
        std::cerr << ", " << secDir << ", secSource=";
        writeLeafPathElementVector(std::cerr, secSource);
        std::cerr << ")" << std::endl;
    }

    if (firSource.empty()) {
        LPECollectionConstIt secIt = secSource.begin();
        for (unsigned ind = 0; secIt != secSource.end(); ++secIt, ++ind) {
            LeafPathElement newLPE(*secIt, secDir, ind);
            peList.push_back(newLPE);
        }
        return;
    }
    if (secSource.empty()) {
        LPECollectionConstIt secIt = firSource.begin();
        for (unsigned ind = 0; secIt != firSource.end(); ++secIt, ++ind) {
            LeafPathElement newLPE(*secIt, firDir, ind);
            peList.push_back(newLPE);
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
        LeafPathElement newLPE(*nLPE, d, nextInd);
        peList.push_back(newLPE);
    }

    if (gDebugging) {
        std::cerr << "mergePathElementLists complete peList =";
        writeLeafPathElementVector(std::cerr, peList);
        std::cerr  << std::endl;
    }

}


inline long CentroidScore(const NxsSimpleNode &nd,  long numActive) {
    NdBlob * b = CurrentBlob(nd);
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
        EdgeDecompInfo * edgeInfo = (CurrentBlob(nd))->GetActiveEdgeInfoPtr();
        assert(edgeInfo);
        if (dir == LEFT_OR_RIGHT) {
            LPECollection &peList = edgeInfo->GetClosestLeavesAboveRef();
            assert(peList.empty());
            mergePathElementLists(peList, LEFT_DIR, GetClosestLeavesAboveNd(&o),
                                          RIGHT_DIR, GetClosestLeavesAboveNd(&o2));
            edgeInfo->SetAboveInitialized(true);
        }
        else {
            const LPECollection & parSource = GetClosestLeavesBelowNd(&o2);
            LPECollection &peList = edgeInfo->GetClosestLeavesBelowRef();
            TreeDirection od = (dir == LEFT_OR_BELOW ? LEFT_DIR : RIGHT_DIR);
            mergePathElementLists(peList, od, GetClosestLeavesAboveNd(&o), BELOW_DIR, parSource);
            edgeInfo->SetBelowInitialized(true);
        }
    }
}

NdBlobSettingStruct stageSingleNodeFlagActivePathDown(const NxsSimpleNode * currAnc,
                                  const LPECollection * newActiveLeaves,
                                  const NxsSimpleNode * & nextNd,
                                  const LPECollection * & nextLPEC,
                                  const NxsSimpleNode * & futureNd,
                                  const LPECollection * & futureLPEC,
                                  TreeSweepDirection * futureBitToAdd,
                                  TreeSweepDirection bitToAdd,
                                  TreeDirection focalEdgeDirection,
                                  long numActiveAbove,
                                  const NxsSimpleNode * deepestAnc,
                                  std::set<NdBlob *> & alteredBlobs) {
    if (gDebugging) {
        std::cerr << "stageSingleNodeFlagActivePathDown(currAnc=" << (long)currAnc->GetTaxonIndex(); 
        if (currAnc->GetFirstChild()) {
            std::cerr << "LeftChild(*currAnc)=" << (long)currAnc->GetFirstChild()->GetTaxonIndex();
            if (RightChild(*currAnc)) {
                std::cerr << ", RightChild(*currAnc)=" << (long)RightChild(*currAnc)->GetTaxonIndex();
            }
            if (currAnc->GetParent()) {
                std::cerr << ", Parent(*currAnc)=" << (long)Parent(*currAnc)->GetTaxonIndex();
            }
        }
        std::cerr << ", numActiveAbove=" << (long) numActiveAbove << ')' << std::endl;
        if (newActiveLeaves) {
            std::cerr << "stageSingleNodeFlagActivePathDown leafset :";
            writeLeafSet(std::cerr, *newActiveLeaves);
            std::cerr << '\n';
        }
    }
    assert(currAnc);
    NdBlob * currBlob = CurrentBlob(*currAnc);
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
    if (currAnc == gRoot || currAnc == deepestAnc) {
        return NdBlobSettingStruct(0L, BELOW_DIR, numActiveAbove, BELOW_DIR_BIT, 0L, alteredBlobs);
    }

    if (false && IsRootChild(*currAnc)) {
        return NdBlobSettingStruct(currBlob, focalEdgeDirection, numActiveAbove, TreeSweepDirection(BELOW_DIR_BIT | bitToAdd), currBlob->GetFullEdgeInfoPtr(), alteredBlobs);
    }
    
    if (newActiveLeaves == 0L) {
        if (currBlob->IsParentsLeftChild()) {
            *futureBitToAdd = RIGHT_DIR_BIT;
            nextNd = RightChild(*Parent(*currAnc));
        }
        else {
            *futureBitToAdd = LEFT_DIR_BIT;
            nextNd = LeftChild(*Parent(*currAnc));
        }
        nextLPEC = 0L;
        futureNd = Parent(*currAnc);
        futureLPEC = 0L;
        return NdBlobSettingStruct(currBlob, focalEdgeDirection, numActiveAbove, TreeSweepDirection(BELOW_DIR_BIT | bitToAdd), currBlob->GetFullEdgeInfoPtr(), alteredBlobs);
    }

    LPECollectionConstIt nalIt = newActiveLeaves->begin();
    for (; nalIt != newActiveLeaves->end(); ++nalIt) {
        int indInCurr = nalIt->indexInNext;
        const LeafPathElement & el = currBlob->GetFullEdgeInfoConstPtr()->GetCloseLeavesBelow().at(indInCurr);
        if (gDebugging) {
            std::cerr << "nalIt->indexInNext = " << indInCurr << " el.leaf=" << el.GetLeaf()->GetTaxonIndex() << " nalIt->leaf=" << (long) nalIt->GetLeaf()->GetTaxonIndex() << " el.dirToNext = " << el.dirToNext  << std::endl;
            std::cerr << "newActiveLeaves: ";
            writeLeafSet(std::cerr, *newActiveLeaves);
            std::cerr << '\n';
            std::cerr << "currBlob->GetFullEdgeInfoConstPtr()->GetCloseLeavesBelow(): ";
            writeLeafSet(std::cerr, currBlob->GetFullEdgeInfoConstPtr()->GetCloseLeavesBelow());
            std::cerr << '\n';
        }
        assert(el.GetLeaf() == nalIt->GetLeaf());
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

    if (currBlob->lpeScratch1.empty()) {
        if (currBlob->lpeScratch2.empty()) {
             //INIT?
            return NdBlobSettingStruct(currBlob, focalEdgeDirection, numActiveAbove, bitToAdd, new EdgeDecompInfo(__FILE__, __LINE__), alteredBlobs);
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
    EdgeDecompInfo * edi = new EdgeDecompInfo(__FILE__, __LINE__);  //INIT?
    
    return NdBlobSettingStruct(currBlob, focalEdgeDirection, numActiveAbove, TreeSweepDirection(BELOW_DIR_BIT | bitToAdd), edi, alteredBlobs);
}


void NdBlobSettingStruct::SetObject() {
    assert(blob);
    if (blob == 0L)
        return;
    if (gDebugging) {
        std::cerr << "SetObject before ";
        debugBlobPrint(std::cerr, blob);
        std::cerr << std::endl;
    }
    assert(this->alteredBlobs.find(blob) == this->alteredBlobs.end());
    this->alteredBlobs.insert(blob);
    blob->SetActiveLeafDir(activeLeafDir);
    blob->SetFocalEdgeDir(focalDir);
    blob->SetNumActiveLeavesAbove(numActiveAbove);
    blob->SetActiveEdgeInfoPtr(edgeDecompInfo);

    if (gDebugging) {
        std::cerr << "SetObject after ";
        debugBlobPrint(std::cerr, blob);
        std::cerr << std::endl;
    }

    blob = 0L; // we don't want to accidentally push two sets of variables onto the stack
}

void nonRecursiveFlagActivePathDown(const NxsSimpleNode * currAnc,
                                  const LPECollection * newActiveLeaves,
                                  const NxsSimpleNode * & nextNd,
                                  const LPECollection * & nextLPEC,
                                  const NxsSimpleNode * & futureNd,
                                  const LPECollection * & futureLPEC,
                                  TreeSweepDirection * futureBitToAdd,
                                  TreeSweepDirection bitToAdd,
                                  TreeDirection focalEdgeDirection,
                                  long numActiveAbove,
                                  const NxsSimpleNode * deepestAnc,
                                  std::set<NdBlob *> & alteredBlobs) {
    NdBlobSettingStruct bss = stageSingleNodeFlagActivePathDown(currAnc, newActiveLeaves,
                                  nextNd,
                                  nextLPEC,
                                  futureNd,
                                  futureLPEC,
                                  futureBitToAdd,
                                  bitToAdd,
                                  focalEdgeDirection, 
                                  numActiveAbove,
                                  deepestAnc,
                                  alteredBlobs);

    if (bss.blob) {
        bss.SetObject();
    }
}
NdBlobSettingStruct stageNonRecursiveFlagActivePathUp(const NxsSimpleNode * currAnc,
                                  const LPECollection * newActiveLeaves,
                                  const NxsSimpleNode * & nextNd,
                                  const LPECollection * & nextLPEC,
                                  const NxsSimpleNode * & futureNd,
                                  const LPECollection * & futureLPEC,
                                  std::set<NdBlob *> & alteredBlobs) {
    if (gDebugging) {
        std::cerr << "nonRecursiveFlagActivePathUp(currAnc=" << (long)currAnc->GetTaxonIndex() ;
        if (currAnc->GetFirstChild()) {
            std::cerr << "LeftChild(*currAnc)=" << (long)currAnc->GetFirstChild()->GetTaxonIndex();
            std::cerr << ", RightChild(*currAnc)=" << (long)RightChild(*currAnc)->GetTaxonIndex();
        }
        std::cerr << ", newActiveLeaves=" << (long) newActiveLeaves << ')' << std::endl;
    }
    assert(currAnc);
    NdBlob * currBlob = CurrentBlob(*currAnc);
    assert(currBlob);

    currBlob->lpeScratch1.clear();
    currBlob->lpeScratch2.clear();
    nextNd = 0L;
    nextLPEC = 0L;
    futureNd = 0L;
    futureLPEC = 0L;
    if (currAnc->IsTip() || newActiveLeaves == 0L) {
        return NdBlobSettingStruct(0L, BELOW_DIR, 0L, BELOW_DIR_BIT, 0L, alteredBlobs);
    }

    LPECollectionConstIt nalIt = newActiveLeaves->begin();
    LPECollection closeAbove;
    for (; nalIt != newActiveLeaves->end(); ++nalIt) {
        int indInCurr = nalIt->indexInNext;
        const LeafPathElement & el = currBlob->GetFullEdgeInfoConstPtr()->GetCloseLeavesAbove().at(indInCurr);
        if (gDebugging) {
            std::cerr << "nalIt->indexInNext = " << indInCurr << " el.leaf=" << el.GetLeaf()->GetTaxonIndex() << " nalIt->leaf=" << nalIt->GetLeaf()->GetTaxonIndex() << std::endl;
        }
        assert(el.GetLeaf() == nalIt->GetLeaf());
        if (el.dirToNext == LEFT_DIR) {
            currBlob->lpeScratch1.push_back(el);
        }
        else {
            assert(el.dirToNext == RIGHT_DIR);
            currBlob->lpeScratch2.push_back(el);
        }
        closeAbove.push_back(el);
    }
    
    TreeSweepDirection tsd = NO_DIR_BIT;
    if (currBlob->lpeScratch1.empty()) {
        if (currBlob->lpeScratch2.empty()) {
            return NdBlobSettingStruct(0L, BELOW_DIR, 0L, BELOW_DIR_BIT, 0L, alteredBlobs);
        }
        tsd = RIGHT_DIR_BIT;
        futureNd = RightChild(*currAnc);
        futureLPEC = &(currBlob->lpeScratch2);
    }
    else if (currBlob->lpeScratch2.empty()) {
        tsd = LEFT_DIR_BIT;
        nextNd = LeftChild(*currAnc);
        nextLPEC = &(currBlob->lpeScratch1);
        }
    else {
        tsd = LEFT_OR_RIGHT;
        nextNd = LeftChild(*currAnc);
        nextLPEC = &(currBlob->lpeScratch1);
        futureNd = RightChild(*currAnc);
        futureLPEC = &(currBlob->lpeScratch2);
    }
    
    EdgeDecompInfo * edi = new EdgeDecompInfo(__FILE__, __LINE__); 
    edi->SetClosestLeavesAbove(closeAbove);
    return NdBlobSettingStruct(currBlob, BELOW_DIR, newActiveLeaves->size(), tsd, edi, alteredBlobs);
}



void nonRecursiveFlagActivePathUp(const NxsSimpleNode * currAnc,
                                  const LPECollection * newActiveLeaves,
                                  const NxsSimpleNode * & nextNd,
                                  const LPECollection * & nextLPEC,
                                  const NxsSimpleNode * & futureNd,
                                  const LPECollection * & futureLPEC,
                                  std::set<NdBlob*> &alteredBlobs) {
    NdBlobSettingStruct bss =  stageNonRecursiveFlagActivePathUp(currAnc,
                                  newActiveLeaves,
                                  nextNd,
                                  nextLPEC,
                                  futureNd,
                                  futureLPEC,
                                  alteredBlobs);
    if (bss.blob) {
        bss.SetObject();
    }
}
// Move up the tree.. updating the NdBlob for currAnc and all of its descendants.
// Pushing current activeEdgeInfo and activeLeafDir to their stacks and update based
//  on the `newActiveLeaves`  If `newActiveLeaves` is 0L, then the entire subtree
//  should be regarded as active
/// Sets focalEdgeDir for nodes along the path *EXCEPT* for currAnc
void flagActivePathUp(const NxsSimpleNode *currAnc, const LPECollection *newActiveLeaves, std::set<NdBlob*> &alteredBlobs) {
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
        nonRecursiveFlagActivePathUp(currNd, currLPEC, nextNd, nextLPEC, futureNd, futureLPEC, alteredBlobs);

        if (futureNd != 0) {
            ndStack.push(futureNd);
            lpecStack.push(futureLPEC);

            NdBlob * futureBlob = CurrentBlob(*futureNd);
            assert(futureBlob);
        }
        if (nextNd != 0) {
            ndStack.push(nextNd);
            lpecStack.push(nextLPEC);

            NdBlob * nextBlob = CurrentBlob(*nextNd);
            assert(nextBlob);
        }
    }
}


/// Sets focalEdgeDir for nodes along the path *EXCEPT* for currAnc
/// \returns the "active Root"
const NxsSimpleNode * flagActivePathDown(
  const NxsSimpleNode *currAnc,
  const LPECollection *newActiveLeaves,
  TreeSweepDirection dirOfOtherActive, 
  long numAboveDeepest,
  const NxsSimpleNode * deepestAnc,
  std::set<NdBlob*> &alteredBlobs) {
    if (gDebugging) {
        std::cerr << "flagActivePathDown(currAnc=" << (long) currAnc << ", newActiveLeaves=" <<(long)newActiveLeaves << ", numAboveDeepest=" <<  numAboveDeepest << ", deepestAnc=" << deepestAnc->GetTaxonIndex() << ")" << std::endl;
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
    std::stack<TreeDirection> downFocalDirStack;
    const NxsSimpleNode * actRoot = currAnc;
    TreeDirection actRootFocalDir = LEFT_DIR;
    TreeSweepDirection futureBitToAdd = dirOfOtherActive;
    downStack.push(currAnc);
    downLPECStack.push(newActiveLeaves);
    downBitStack.push(futureBitToAdd);
    assert(dirOfOtherActive == LEFT_DIR_BIT || dirOfOtherActive == RIGHT_DIR_BIT);
    downFocalDirStack.push(dirOfOtherActive == LEFT_DIR_BIT ? RIGHT_DIR : LEFT_DIR);

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
            const TreeDirection currDirToFocal = downFocalDirStack.top();
            downFocalDirStack.pop();
            const NxsSimpleNode * nextNd = 0L;
            const NxsSimpleNode * futureNd = 0L;
            const LPECollection * nextLPEC = 0L;
            const LPECollection * futureLPEC = 0L;
            nonRecursiveFlagActivePathDown(currNd, currLPEC, nextNd, nextLPEC, futureNd, futureLPEC, &futureBitToAdd, currBitToAdd, currDirToFocal, numAboveDeepest, deepestAnc, alteredBlobs);
            NdBlob * currBlob = CurrentBlob(*currNd);
            
            if (futureNd != 0) {
                downStack.push(futureNd);
                downLPECStack.push(futureLPEC);
                downBitStack.push(futureBitToAdd);
                downFocalDirStack.push((currBlob->IsParentsLeftChild() ? LEFT_DIR : RIGHT_DIR));
            }
            if (nextNd != 0) {
                actRoot = Parent(*nextNd);
                assert(Parent(*currNd) == actRoot);
                actRootFocalDir = (currBlob->IsParentsLeftChild() ? LEFT_DIR : RIGHT_DIR);
                
                upStack.push(nextNd);
                upLPECStack.push(nextLPEC);
                if (nextLPEC == 0L) {
                    numAboveDeepest += CurrentBlob(*nextNd)->GetNumActiveLeavesAbove();
                }
                else {
                    numAboveDeepest += nextLPEC->size();
                }
                
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
            nonRecursiveFlagActivePathUp(currNd, currLPEC, nextNd, nextLPEC, futureNd, futureLPEC, alteredBlobs);
            
            if (futureNd != 0) {
                upStack.push(futureNd);
                upLPECStack.push(futureLPEC);
            }
            if (nextNd != 0) {
                upStack.push(nextNd);
                upLPECStack.push(nextLPEC);
            }
        }
    }
    NdBlob * actBlob = CurrentBlob(*actRoot);
    assert(actBlob);
    assert(!actRoot->IsTip());
    NdBlob * leftBlob = CurrentBlob(*LeftChild(*actRoot));
    NdBlob * rightBlob = CurrentBlob(*RightChild(*actRoot)) ;


    EdgeDecompInfo * edi = new EdgeDecompInfo(__FILE__, __LINE__);
    mergePathElementLists(edi->GetClosestLeavesAboveRef(), LEFT_DIR, GetClosestLeavesAboveNd(LeftChild(*actRoot)),
                                          RIGHT_DIR, GetClosestLeavesAboveNd(RightChild(*actRoot)));
    edi->SetAboveInitialized(true);
    NdBlobSettingStruct nbss(actBlob, 
                             actRootFocalDir, 
                             leftBlob->GetNumActiveLeavesAbove() + rightBlob->GetNumActiveLeavesAbove(), 
                             LEFT_OR_RIGHT, 
                             edi, 
                             alteredBlobs);
    nbss.SetObject();
    
    return actRoot;
}

void writeNewickOfActive(std::ostream & out, const NxsSimpleNode *activeRoot) {
    assert(activeRoot);
    typedef std::pair<const NxsSimpleNode *, std::stack<char> > NdSuffixPair;
    std::stack<NdSuffixPair> ndStack;
    std::stack<char> currSuffix;
    char currPrefix = '\0';
    for (;;) {
        NdBlob * currBlob = CurrentBlob(*activeRoot);
        //if (gDebugging)
        //    std::cerr << "activeRoot = " << (long) activeRoot << ", activeRoot->GetTaxonIndex()=" << activeRoot->GetTaxonIndex() << ", currBlob->activeLeafDir="  << currBlob->activeLeafDir << std::endl;
        if (activeRoot->IsTip()) {
            out << (activeRoot->GetTaxonIndex() + 1);
            activeRoot = 0L;
        }
        else {
            TreeSweepDirection fullActiveLeafDir = (currBlob ? currBlob->GetActiveAndFocalLeafDir() : LEFT_OR_RIGHT);
            //if (gDebugging)
            //    std::cerr << " fullActiveLeafDir = " << fullActiveLeafDir << std::endl;
            const bool leftActive = fullActiveLeafDir & LEFT_DIR_BIT;
            const bool rightActive = fullActiveLeafDir & RIGHT_DIR_BIT;
            if (leftActive && rightActive) {
                if (currPrefix != '\0') {
                    currSuffix.push(currPrefix);
                }
#               if defined(DEBUGGING_OUTPUT)
                    if (gDebugging) {
                        out << "[i" << activeRoot->GetTaxonIndex() << ']';
                    }
#               endif
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

const NxsSimpleNode * findCentroidForActiveLeavesUp(const NxsSimpleNode *nd,
                                                  long totalNumActiveLeaves,
                                                  long ndScore,
                                                  long * retScore) {
    if (nd->IsTip()) {
        if (gDebugging)
            std::cerr << "findCentroidForActiveLeavesUp: Returning 0L from tip\n";
        return 0L;
    }
    if (gDebugging)
        std::cerr << "findCentroidForActiveLeavesUp: ndScore = " << ndScore << "\n";
    
    NdBlob * currBlob = CurrentBlob(*nd);
    *retScore = ndScore;
    const NxsSimpleNode * possibleReturn;
    TreeSweepDirection fullActiveLeafDir = currBlob->GetActiveAndFocalLeafDir();
    if (!nd->IsTip()) {
        const NxsSimpleNode * child;
        long childScore;
        if (fullActiveLeafDir & LEFT_DIR_BIT) {
            if (gDebugging)
                std::cerr << "findCentroidForActiveLeavesUp: Calculating LeftChild centroid score... ";
            child = LeftChild(*nd);
            childScore = CentroidScore(*child, totalNumActiveLeaves);
            if (childScore <= ndScore) {
                if (gDebugging)
                    std::cerr << "findCentroidForActiveLeavesUp: Recursing left from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << childScore << '\n';
                possibleReturn = findCentroidForActiveLeavesUp(child, totalNumActiveLeaves, childScore, retScore);
                if (possibleReturn) {
                    if (*retScore < ndScore) {
                        return possibleReturn;
                    }
                }
                else if (childScore < ndScore) {
                    *retScore = childScore;
                    if (gDebugging)
                        std::cerr << "findCentroidForActiveLeavesUp: Returning Left from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << *retScore << '\n';
                    return child;
                }
            }
        }
        if (fullActiveLeafDir & RIGHT_DIR_BIT) {
            if (gDebugging)
                std::cerr << "findCentroidForActiveLeavesUp: Calculating RightChild centroid score... ";
            child = RightChild(*nd);
            childScore = CentroidScore(*child, totalNumActiveLeaves);
            if (childScore <= ndScore) {
                if (gDebugging)
                    std::cerr << "findCentroidForActiveLeavesUp: Recursing right from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << childScore << '\n';
                possibleReturn = findCentroidForActiveLeavesUp(child, totalNumActiveLeaves, childScore, retScore);
                if (possibleReturn) {
                    if (*retScore < ndScore) {
                        return possibleReturn;
                    }
                }
                else if (childScore < ndScore) {
                    *retScore = childScore;
                    if (gDebugging)
                        std::cerr << "findCentroidForActiveLeavesUp: Returning Right from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << *retScore << '\n';
                    return child;
                }
            }
        }
    }
    if (gDebugging)
        std::cerr << "findCentroidForActiveLeavesUp: Centroid Edge found nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ")\n";
    return nd;
}

///
// Moves through down through ancestors of des until finding a node for which
//      GetActiveAndFocalLeafDir has both the LEFT_DIR_BIT and RIGHT_DIR_BIT 
//      set.
//  returns the pair that is the ancestor, and the child of the ancestor that is
//      not along the path from the ancestor to `des`.
///
std::pair<const NxsSimpleNode *, const NxsSimpleNode *> NextAncAndOtherDesWithBothChildrenActive(const NxsSimpleNode &des) {
    const NxsSimpleNode * prev = &des;
    const NxsSimpleNode * nd = Parent(des);
    for (;;) {
        NdBlob * currBlob = CurrentBlob(*nd);
        TreeSweepDirection fullActiveLeafDir = currBlob->GetActiveAndFocalLeafDir();
        if (TreeSweepDirection(fullActiveLeafDir & LEFT_OR_RIGHT) == LEFT_OR_RIGHT) {
            const NxsSimpleNode * other = ((CurrentBlob(*prev))->IsParentsLeftChild() ? RightChild(*nd) : LeftChild(*nd));
            return std::pair<const NxsSimpleNode *, const NxsSimpleNode *>(nd, other);
        }
        else if (fullActiveLeafDir & BELOW_DIR_BIT) {
            prev = nd;
            nd = Parent(*nd);
        }
        else
            return std::pair<const NxsSimpleNode *, const NxsSimpleNode *>(0L, 0L);
    }
}
const NxsSimpleNode * findCentroidForActiveLeavesDown(const NxsSimpleNode *nd,
                                                  long totalNumActiveLeaves,
                                                  long ndScore,
                                                  long * retScore) {
    assert(!nd->IsTip());
    if (gDebugging)
        std::cerr << "findCentroidForActiveLeavesDown: ndScore = " << ndScore << "\n";
    
    NdBlob * currBlob = CurrentBlob(*nd);
    *retScore = ndScore;
    const NxsSimpleNode * possibleReturn;
    TreeSweepDirection fullActiveLeafDir = currBlob->GetActiveAndFocalLeafDir();
    if (fullActiveLeafDir & BELOW_DIR_BIT ) {
        std::pair<const NxsSimpleNode *, const NxsSimpleNode *> p = NextAncAndOtherDesWithBothChildrenActive(*nd);
        const NxsSimpleNode * anc = p.first;
        const NxsSimpleNode * other = p.second;
        assert(anc);
        if (gDebugging)
            std::cerr << "findCentroidForActiveLeavesDown: Calculating Anc centroid score... ";
        long ancScore = CentroidScore(*anc, totalNumActiveLeaves);
        long otherScore = CentroidScore(*other, totalNumActiveLeaves);
        if (ancScore <= ndScore) {
            if (otherScore <= ancScore) {
                if (gDebugging)
                    std::cerr << "findCentroidForActiveLeavesDown: Moving to \"sib\" subtree from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << otherScore << '\n';
                possibleReturn = findCentroidForActiveLeavesUp(other, totalNumActiveLeaves, otherScore, retScore); // we move up the sib subtree...
                if (possibleReturn) {
                    if (*retScore < ndScore) {
                        return possibleReturn;
                    }
                }
            }
            else {
                if (gDebugging)
                    std::cerr << "findCentroidForActiveLeavesDown: Recursing down from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << ancScore << '\n';
                possibleReturn =  findCentroidForActiveLeavesDown(nd, totalNumActiveLeaves, ancScore, retScore);
                if (possibleReturn) {
                    if (*retScore < ndScore) {
                        return possibleReturn;
                    }
                }
                else {
                    if (otherScore <= ancScore) {
                        *retScore = otherScore;
                        if (gDebugging)
                            std::cerr << "findCentroidForActiveLeavesDown: Returning Anc from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << *retScore << '\n';
                        return other;
                    }
                    else {
                        *retScore = ancScore;
                        if (gDebugging)
                            std::cerr << "findCentroidForActiveLeavesDown: Returning Anc from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << *retScore << '\n';
                        return anc;
                    }
                }
            }
        }
        else if (otherScore <= ndScore) {
            if (gDebugging)
                std::cerr << "findCentroidForActiveLeavesDown: Moving to \"sib\" subtree from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << otherScore << '\n';
            possibleReturn = findCentroidForActiveLeavesUp(other, totalNumActiveLeaves, otherScore, retScore); // we move up the sib subtree...
            if (possibleReturn) {
                if (*retScore < ndScore) {
                    return possibleReturn;
                }
            }
            else if (otherScore < ndScore) {
                *retScore = otherScore;
                if (gDebugging)
                    std::cerr << "findCentroidForActiveLeavesDown: Returning \"sib\" from nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ") toward score = " << *retScore << '\n';
                return other;
            }
        }
    }
    if (gDebugging)
        std::cerr << "findCentroidForActiveLeavesDown: Centroid Edge found nd[i" << nd->GetTaxonIndex() << "] (score=" << ndScore << ")\n";
    return nd;
}

const NxsSimpleNode * doDecomposition(const NxsSimpleNode * topCentroidChild,
                         const NxsSimpleNode * leftSubtreeRoot, const LPECollection * leftCommonLeafSet,
                         const NxsSimpleNode * rightSubtreeRoot, const LPECollection * rightCommonLeafSet,
                         const NxsSimpleNode * sibSubtreeRoot, const LPECollection * sibCommonLeafSet,
                         const NxsSimpleNode * parSubtreeRoot, const LPECollection * parCommonLeafSet,
                        NxsString namePrefix,
                        long numActiveLeavesPrev,
                        const NxsSimpleNode * effectiveRoot,
                        bool movingUp) {

    NdBlob * centroidBlob = CurrentBlob(*topCentroidChild);
    NdBlob * parBlob = CurrentBlob(*parSubtreeRoot); 
    NdBlob * sibBlob = CurrentBlob(*sibSubtreeRoot); 
    NdBlob * leftBlob = CurrentBlob(*leftSubtreeRoot);
    NdBlob * rightBlob = CurrentBlob(*rightSubtreeRoot);

    const unsigned numLeftChosen = (leftCommonLeafSet ? leftCommonLeafSet->size() : leftBlob->GetNumActiveLeavesAbove());
    const unsigned numSibChosen =  (sibCommonLeafSet ? sibCommonLeafSet->size() : sibBlob->GetNumActiveLeavesAbove());
    const unsigned numBelowChosen = (parCommonLeafSet ? parCommonLeafSet->size() : numActiveLeavesPrev - parBlob->GetNumActiveLeavesAbove());
    const unsigned numRightChosen = (rightCommonLeafSet ? rightCommonLeafSet->size() : rightBlob->GetNumActiveLeavesAbove());
    long sizeOfSubproblem = numBelowChosen + numRightChosen + numLeftChosen + numSibChosen;

    std::set<NdBlob *> alteredBlobs; // \TEMP: by knowing how we do the traversal, we should be able to avoid storing the record of touched nodes in a set...

    if (gDebugging) {
        std::cerr << "Before flagging numLeftChosen = " << numLeftChosen << ", numRightChosen = " << numRightChosen << " numSibChosen = " << numSibChosen << ", numBelowChosen = " << numBelowChosen << std::endl;
        std::cerr << "Common leafset: left: ";
        if (leftCommonLeafSet)
            writeLeafSet(std::cerr, *leftCommonLeafSet);
        else
            std::cerr << "<all leaves>";
        std::cerr << "\nCommon leafset: right: ";
        if (rightCommonLeafSet)
            writeLeafSet(std::cerr, *rightCommonLeafSet);
        else
            std::cerr << "<all leaves>";
        std::cerr << "\nCommon leafset: sib: ";
        if (sibCommonLeafSet)
            writeLeafSet(std::cerr, *sibCommonLeafSet);
        else
            std::cerr << "<all leaves>";
        std::cerr << "\nCommon leafset: parSubtreeRoot: ";
        if (parCommonLeafSet)
            writeLeafSet(std::cerr, *parCommonLeafSet);
        else
            std::cerr << "<all leaves>";
        std::cerr << "\n";
    }

    flagActivePathUp(leftSubtreeRoot, leftCommonLeafSet, alteredBlobs);
    if (gDebugging) {
        std::cerr << "After flagging left\nleftBlob->GetNumActiveLeavesAbove() = " << leftBlob->GetNumActiveLeavesAbove() << ", sizeOfSubproblem = " << sizeOfSubproblem << std::endl;
        std::cerr << "left Node: ";
        debugNdPrint(std::cerr, *leftSubtreeRoot);
        std::cerr << "\n";
        debugPrintTree(std::cerr, leftSubtreeRoot);
        std::cerr << "\n";
    }
    debugCheckActiveFlags(leftSubtreeRoot, numLeftChosen);
    
    flagActivePathUp(rightSubtreeRoot, rightCommonLeafSet, alteredBlobs);
    if (gDebugging) {
        std::cerr << "After flagging right\nright Node: ";
        debugNdPrint(std::cerr, *rightSubtreeRoot);
        std::cerr << "\n";
        debugPrintTree(std::cerr, rightSubtreeRoot);
        std::cerr << "\n";
    }
    debugCheckActiveFlags(rightSubtreeRoot, numRightChosen);

    ////////////////////////////////////////////////////////////////////////////
    // Set the active fields for the centroid path
    EdgeDecompInfo * edi = new EdgeDecompInfo(__FILE__, __LINE__);   //INIT?
    NdBlobSettingStruct centroidNBSS(centroidBlob, BELOW_DIR, numRightChosen + numLeftChosen, LEFT_OR_RIGHT, edi, alteredBlobs);
    centroidNBSS.SetObject();

    const NxsSimpleNode * prev = topCentroidChild;
    const NxsSimpleNode * nd = Parent(*prev);
    while (nd != parSubtreeRoot) {
        if (gDebugging) {
            std::cerr << "Flagging below centroid" << std::endl;
        }
        NdBlob * b = CurrentBlob(*nd);
        NdBlob * prevb = CurrentBlob(*prev);
        TreeSweepDirection d = (prevb->IsParentsLeftChild() ? LEFT_DIR_BIT : RIGHT_DIR_BIT);
        edi = new EdgeDecompInfo(__FILE__, __LINE__);  //INIT?
        edi->SetClosestLeavesAbove(GetClosestLeavesAboveNd(prev)); 
        NdBlobSettingStruct cpBSS(b, BELOW_DIR, numRightChosen + numLeftChosen, d, edi, alteredBlobs);
        cpBSS.SetObject();
        prev = nd;
        nd = Parent(*prev);
    }

    if (gDebugging) {
        std::cerr << "After flagging centroid path." << std::endl;
        std::cerr << "topCentroidChild Node: ";
        debugNdPrint(std::cerr, *rightSubtreeRoot);
        std::cerr << "\n";
        debugPrintTree(std::cerr, SibNode(*sibSubtreeRoot));
        std::cerr << "\n";
    }


    flagActivePathUp(sibSubtreeRoot, sibCommonLeafSet, alteredBlobs);
    if (gDebugging) {
        std::cerr << "After flagging sib\nsib Node: ";
        debugNdPrint(std::cerr, *sibSubtreeRoot);
        std::cerr << "\n";
        debugPrintTree(std::cerr, sibSubtreeRoot);
        std::cerr << "\n";
    }
    debugCheckActiveFlags(sibSubtreeRoot, numSibChosen);


    TreeSweepDirection parOtherActiveDir = sibBlob->IsParentsLeftChild() ? LEFT_DIR_BIT : RIGHT_DIR_BIT;
    long numActiveAboveCentroid = numRightChosen + numLeftChosen;
    long numAboveParInLeftSubproblem = numActiveAboveCentroid + numSibChosen;


    const NxsSimpleNode * activeRoot = flagActivePathDown(parSubtreeRoot, parCommonLeafSet, parOtherActiveDir, numAboveParInLeftSubproblem, effectiveRoot, alteredBlobs);
    if (gDebugging) {
        std::cerr << "After flagging par\npar Node: sizeOfSubproblem = " << sizeOfSubproblem << '\n';
        std::cerr << "parSubtreeRoot => ";
        debugNdPrint(std::cerr, *parSubtreeRoot);
        std::cerr << "\n";
        std::cerr << "activeRoot => ";
        debugNdPrint(std::cerr, *activeRoot);
        std::cerr << "\n";
        std::cerr << "effectiveRoot => ";
        debugNdPrint(std::cerr, *effectiveRoot);
        std::cerr << "\n";
        debugPrintTree(std::cerr, effectiveRoot);
        debugCheckActiveFlags(activeRoot, sizeOfSubproblem);
    }
    

    
        
    
    if (gOutputStream != 0L && (sizeOfSubproblem <= gMaxSubProblemSize || gEmitIntermediates)) {
        *gOutputStream << "Tree " << namePrefix << " = [&U] ";
        writeNewickOfActive(*gOutputStream, activeRoot);
        *gOutputStream << ";\n";
    }
    
    
    
    
    
    if (sizeOfSubproblem > gMaxSubProblemSize) {
        // we need to "fix" the EdgeDecompInfo structs along the centroid path, because they might be used...
        if (movingUp) {
            // we are going to be moving "up" in this subtree. So we'll need the centroid nodes paths to be have correct closestLeavesBelow fields...
            
        }
        else {
            // we are going to be moving "down" in this subtree. So we'll need the centroid nodes paths to be have correct closestLeavesAbove fields...
            const NxsSimpleNode * prev = topCentroidChild;
            
            EdgeDecompInfo * cedi = (CurrentBlob(*prev))->GetActiveEdgeInfoPtr();
            mergePathElementLists(cedi->GetClosestLeavesAboveRef(), LEFT_DIR, GetClosestLeavesAboveNd(leftSubtreeRoot),
                                          RIGHT_DIR, GetClosestLeavesAboveNd(rightSubtreeRoot));
            cedi->SetAboveInitialized(true);
            const NxsSimpleNode * nd = Parent(*prev);
            const LPECollection emptyCollection;
            while (nd != parSubtreeRoot) {
                cedi = CurrentBlob(*nd)->GetActiveEdgeInfoPtr();
                NdBlob * prevb = CurrentBlob(*prev);
                if (prevb->IsParentsLeftChild()) {
                    mergePathElementLists(cedi->GetClosestLeavesAboveRef(), LEFT_DIR, GetClosestLeavesAboveNd(prev),
                                          RIGHT_DIR, emptyCollection);
                }
                else {
                    mergePathElementLists(cedi->GetClosestLeavesAboveRef(), LEFT_DIR, emptyCollection,
                                          RIGHT_DIR, GetClosestLeavesAboveNd(prev));
                }
                cedi->SetAboveInitialized(true);
                nd = Parent(*prev);
            }
        }
    

    
    
        long sc = CentroidScore(*rootLeft, sizeOfSubproblem);
        long minCS;
        const NxsSimpleNode * nextCentroid = 0L;
        if (leftCommonLeafSet == 0L || rightCommonLeafSet == 0L) {
            nextCentroid = findCentroidForActiveLeavesUp(topCentroidChild, sizeOfSubproblem, sc, &minCS);
        }
        else if (sibCommonLeafSet == 0L || parCommonLeafSet == 0L) {
            nextCentroid = findCentroidForActiveLeavesDown(topCentroidChild, sizeOfSubproblem, sc, &minCS);
        }
        else {
            assert(false); // in the current system of decompositions this should not occur.
            nextCentroid = findCentroidForActiveLeavesUp(activeRoot, sizeOfSubproblem, sc, &minCS);
        }
        if (gDebugging) {
            std::cerr << "nextCentroid: " ; debugNdPrint(std::cerr, *nextCentroid); std::cerr << std::endl;
        }
        
        assert(nextCentroid);
        const NxsSimpleNode * nextBottomCentroid = nextCentroid;
        for (;;) {
            const NxsSimpleNode * parnb  = Parent(*nextBottomCentroid);
            NdBlob * currBlob = CurrentBlob(*parnb);
            TreeSweepDirection fullActiveLeafDir = currBlob->GetActiveAndFocalLeafDir();
            if ((fullActiveLeafDir & LEFT_DIR_BIT) && (fullActiveLeafDir & RIGHT_DIR_BIT)) {
                break;
            }
            nextBottomCentroid = parnb;
        }

        const NxsSimpleNode * nextTopCentroid = nextCentroid;
        for (;;) {
            NdBlob * currBlob = CurrentBlob(*nextTopCentroid);
            TreeSweepDirection fullActiveLeafDir = currBlob->GetActiveLeafDir();
            if (fullActiveLeafDir == LEFT_DIR_BIT) {
                nextTopCentroid = LeftChild(*nextTopCentroid);
            }
            else if (fullActiveLeafDir == RIGHT_DIR_BIT) {
                nextTopCentroid = RightChild(*nextTopCentroid);
            }
            else if (fullActiveLeafDir == BELOW_DIR_BIT) {
                nextTopCentroid = Parent(*nextTopCentroid);
            }
            else {
                assert(fullActiveLeafDir == LEFT_OR_RIGHT || fullActiveLeafDir == LEFT_OR_BELOW || fullActiveLeafDir == RIGHT_OR_BELOW || fullActiveLeafDir == ALL_DIR_BITS);
                break;
            }

        }
        decomposeAroundCentroidChild(nextTopCentroid, nextBottomCentroid, sizeOfSubproblem, namePrefix, activeRoot);
    }
    
    
    if (gDebugging) {
        std::cerr << "about to pop attributes from altered blobs..." << std::endl;
    }
    std::set<NdBlob *>::const_iterator abIt = alteredBlobs.begin();
    for (; abIt != alteredBlobs.end(); ++abIt) {
        (*abIt)->popBlob();
    }
    if (gDebugging) {
        std::cerr << "exiting doDecomposition" << std::endl;
    }
    
    return activeRoot;
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
void decomposeAroundCentroidChild(const NxsSimpleNode *topCentroidChild,
                                  const NxsSimpleNode *bottomCentroidChild,
                                  long numActiveLeaves,
                                  NxsString namePrefix,
                                  const NxsSimpleNode * effectiveRoot) {
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
    const NdBlob * centroidBlob = CurrentBlob(*topCentroidChild);
    const NdBlob * parBlob = CurrentBlob(*parSubtreeRoot); 
    const NdBlob * sibBlob = CurrentBlob(*sibSubtreeRoot); 
    const NdBlob * leftBlob = CurrentBlob(*leftSubtreeRoot);
    const NdBlob * rightBlob = CurrentBlob(*rightSubtreeRoot);
    
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

    const long floorHalfAbove = numAboveChosen/2;
    const long cielHalfAbove = (2*floorHalfAbove == numAboveChosen ? floorHalfAbove : 1 + floorHalfAbove);
    assert(floorHalfAbove + cielHalfAbove == numAboveChosen);
    const long origNumActiveLeavesAboveLeft = leftBlob->GetNumActiveLeavesAbove();
    if (origNumActiveLeavesAboveLeft < floorHalfAbove) {
        numLeftChosen = origNumActiveLeavesAboveLeft;
        assert(rightBlob->GetNumActiveLeavesAbove() >= floorHalfAbove - numLeftChosen);
        numRightChosen = numAboveChosen - numLeftChosen;
    }
    else {
        if (rightBlob->GetNumActiveLeavesAbove() < cielHalfAbove) {
            numRightChosen = rightBlob->GetNumActiveLeavesAbove();
            assert(origNumActiveLeavesAboveLeft >= numAboveChosen - numRightChosen);

        }
        else {
            numRightChosen = cielHalfAbove; //\TODO should probably sort the PLE's and use their order to decide which direction gets the rounding error...
        }
        numLeftChosen = numAboveChosen - numRightChosen;
    }
    if (gDebugging)
        std::cerr << "numAboveChosen=" <<numAboveChosen << ", numBelowChosen=" << numBelowChosen << " sibBlob->GetNumActiveLeavesAbove()=" << sibBlob->GetNumActiveLeavesAbove() <<  ", numActiveLeaves=" << numActiveLeaves << ", parBlob->GetNumActiveLeavesAbove()=" << parBlob->GetNumActiveLeavesAbove() << '\n';


    const long floorHalfBelow = numBelowChosen/2;
    const long cielHalfBelow = (floorHalfBelow*2 == numBelowChosen ? floorHalfBelow : 1 + floorHalfBelow);
    assert(floorHalfBelow + cielHalfBelow == numBelowChosen);
    if (sibBlob->GetNumActiveLeavesAbove() < floorHalfBelow) {
        numSibChosen = sibBlob->GetNumActiveLeavesAbove();
        assert((numActiveLeaves  - parBlob->GetNumActiveLeavesAbove()) >= numBelowChosen - numSibChosen);
        numParChosen = numBelowChosen - numSibChosen;
    }
    else {
        if ((numActiveLeaves  - parBlob->GetNumActiveLeavesAbove()) < cielHalfBelow) {
            numParChosen = (numActiveLeaves  - parBlob->GetNumActiveLeavesAbove());
            assert(sibBlob->GetNumActiveLeavesAbove() == numBelowChosen - numParChosen);
        }
        else {
            numParChosen = cielHalfBelow; //\TODO should probably sort the PLE's and use their order to decide which direction gets the rounding error...
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

    if (gDebugging) {
        std::cerr << "floorHalfBelow = " << floorHalfBelow;
        std::cerr << "numLeftChosen = " << numLeftChosen << " leftSubtreeRoot = ";
        debugNdPrint(std::cerr, *leftSubtreeRoot);
        std::cerr << "\nnumRightChosen = " << numRightChosen << " rightSubtreeRoot = ";
        debugNdPrint(std::cerr, *rightSubtreeRoot);
        std::cerr << "\nnumSibChosen = " << numSibChosen << " sibSubtreeRoot = ";
        debugNdPrint(std::cerr, *sibSubtreeRoot);
        std::cerr << "\nnumParChosen = " << numParChosen << " parSubtreeRoot = ";
        debugNdPrint(std::cerr, *parSubtreeRoot);
        std::cerr << std::endl;
    }

    const LPECollection & leftFullLPEC = GetClosestLeavesAboveNd(leftSubtreeRoot);
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
    const LPECollection & rightFullLPEC = GetClosestLeavesAboveNd(rightSubtreeRoot);
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
    const LPECollection & sibFullLPEC = GetClosestLeavesAboveNd(sibSubtreeRoot);
    if (gDebugging) {
        std::cerr << "sibBlob: "; debugBlobPrint(std::cerr, sibBlob); std::cerr << std::endl;
    }    
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
    const LPECollection & parFullLPEC = GetClosestLeavesBelowNd(parSubtreeRoot);
    assert(parFullLPEC.size() >= (unsigned long)numParChosen);
    LPECollection::const_iterator parfLPECIt = parFullLPEC.begin();
    LPECollection parCommonLeafSet;
    parCommonLeafSet.reserve(numParChosen);
    for (unsigned i=0; i < numParChosen; ++i, ++parfLPECIt) {
        parCommonLeafSet.push_back(LeafPathElement(*parfLPECIt, BELOW_DIR, i));
    }
    assert(parCommonLeafSet.size() == (unsigned long)numParChosen);


    if (gDebugging) {
        std::cerr << "numLeftChosen = " << numLeftChosen << ", numRightChosen = " << numRightChosen  << ", numSibChosen = " << numSibChosen << ",  numParChosen = " << numParChosen << std::endl;
    }
    ////////////////////////////////////////////////////////////////////////////
    // Step 3 - 6 activate each subproblem in turn and recurse...
    //
    ////////////////////////////////////////////////////////////////////////////
    if (false && gEmitIntermediates && gOutputStream != 0L) {
        *gOutputStream << "fTree " << namePrefix << " = [&U] ";
        writeNewickOfActive(*gOutputStream, Parent(*rootLeft));
        *gOutputStream << ";\n";
    }
    NxsString n = namePrefix;
    n += ".0";

    if (gDebugging) {
        std::cerr << "Before left decomposition\n";
        debugPrintTree(std::cerr, effectiveRoot);
    }
    doDecomposition(topCentroidChild,
                        leftSubtreeRoot, 0L,
                        rightSubtreeRoot, &rightCommonLeafSet,
                        sibSubtreeRoot, &sibCommonLeafSet,
                        parSubtreeRoot, &parCommonLeafSet,
                        n,
                        numActiveLeaves,
                        effectiveRoot,
                        true);

    if (gDebugging) {
        std::cerr << "after left, before right decomposition\n";
        debugPrintTree(std::cerr, effectiveRoot);
        debugCheckActiveFlags(effectiveRoot, numActiveLeaves);
    }
    n = namePrefix;
    n += ".1";

    doDecomposition(topCentroidChild,
                        leftSubtreeRoot, &leftCommonLeafSet,
                        rightSubtreeRoot, 0L,
                        sibSubtreeRoot, &sibCommonLeafSet,
                        parSubtreeRoot, &parCommonLeafSet,
                        n,
                        numActiveLeaves,
                        effectiveRoot,
                        true);

    if (gDebugging) {
        std::cerr << "after right, before sib decomposition\n";
        debugPrintTree(std::cerr, effectiveRoot);
        debugCheckActiveFlags(effectiveRoot, numActiveLeaves);
    }
    n = namePrefix;
    n += ".2";

    doDecomposition(topCentroidChild,
                        leftSubtreeRoot, &leftCommonLeafSet,
                        rightSubtreeRoot, &rightCommonLeafSet,
                        sibSubtreeRoot, 0L,
                        parSubtreeRoot, &parCommonLeafSet,
                        n,
                        numActiveLeaves,
                        effectiveRoot,
                        false);
    if (gDebugging) {
        std::cerr << "after sib, before parent decomposition\n";
        debugPrintTree(std::cerr, effectiveRoot);
        debugCheckActiveFlags(effectiveRoot, numActiveLeaves);
    }
    n = namePrefix;
    n += ".3";
    doDecomposition(topCentroidChild,
                        leftSubtreeRoot, &leftCommonLeafSet,
                        rightSubtreeRoot, &rightCommonLeafSet,
                        sibSubtreeRoot, &sibCommonLeafSet,
                        parSubtreeRoot, 0L,
                        n,
                        numActiveLeaves,
                        effectiveRoot,
                        false);
    if (gDebugging) {
        std::cerr << "after parent decomposition\n";
        debugPrintTree(std::cerr, effectiveRoot);
        debugCheckActiveFlags(effectiveRoot, numActiveLeaves);
    }
}

void uppassSettingNodesBelowFields(const std::vector<const NxsSimpleNode *> & preorderTraversal) {
    if (gDebugging) {
        std::cerr << "Starting uppass" <<std::endl;
    }

    unsigned currInternalIndex = gNumLeaves + 1; // + 1 for gRoot


    if (CurrentBlob(*gRoot)->numLeavesAboveEdge != gNumLeaves) {
        gErrorMessage << "Expected " << gNumLeaves << " leaves at the root, but only found " << (CurrentBlob(*gRoot))->numLeavesAboveEdge;
        throw NxsException(gErrorMessage);
    }
    
    std::vector<const NxsSimpleNode *>::const_iterator preIt = preorderTraversal.begin();
    
    for (; preIt != preorderTraversal.end(); ++preIt) {
        const NxsSimpleNode * nd = *preIt;
        assert(nd != 0L);
        if (!nd->IsTip()) {
            (const_cast<NxsSimpleNode *>(nd))->SetTaxonIndex(currInternalIndex++);
    
            NdBlob* ndBlob = CurrentBlob(*nd);
            EdgeDecompInfo * edgeInfo = ndBlob->GetActiveEdgeInfoPtr();
            assert(edgeInfo);
            const NxsSimpleNode * par = Parent(*nd);
            if (ndBlob->IsParentsLeftChild())
                mergeShortLeavesLists(*nd, RIGHT_OR_BELOW, *RightChild(*par), *par);
            else
                mergeShortLeavesLists(*nd, LEFT_OR_BELOW, *LeftChild(*par), *par);
            if (IsRootChild(*nd)) {
                assert(edgeInfo->GetCloseLeavesBelow().size() == 1);
                if (gDebugging) {
                    writeLeafPathElementVector(std::cerr, edgeInfo->GetCloseLeavesBelow()); std::cerr << '\n';
                    writeLeafPathElementVector(std::cerr, GetClosestLeavesAboveNd(gRoot)); std::cerr << '\n';
                    writeLeafPathElementVector(std::cerr, GetClosestLeavesAboveNd(LeftChild(*gRoot))); std::cerr << '\n';
                }
                assert(edgeInfo->GetCloseLeavesBelow()[0].GetScore() == 1);
            }
        }
        if (gDebugging) {
            NdBlob * nb = CurrentBlob(*nd);
            std::cerr << "nd = " << (long) nd << " blob=" << (long) nb << " blob->numLeavesAboveEdge = " << nb->numLeavesAboveEdge <<  " blob->GetActiveEdgeInfoConstPtr() = " << (long) nb->GetActiveEdgeInfoConstPtr() <<  " blob->GetActiveEdgeInfoConstPtr()->GetClosestLeavesAbove().size() = " << GetClosestLeavesAboveNd(nd).size() << '\n';
            if (nd->IsTip()) {
                std::cerr << "nd is tip, so below info not calculated" << std::endl;
            }
            else {
                std::cerr << "nd = " << (long) nd << " blob=" << (long) nb << " blob->numLeavesBelowEdge = " << (gNumLeaves - nb->numLeavesAboveEdge) <<  " blob->GetActiveEdgeInfoConstPtr() = " << (long) nb->GetActiveEdgeInfoConstPtr() <<  " blob->GetActiveEdgeInfoConstPtr()->GetCloseLeavesBelow().size() = " << nb->GetActiveEdgeInfoConstPtr()->GetCloseLeavesBelow().size() << '\n';
            }
        }
    }
    if (gDebugging) {
        std::cerr << "Starting uppass" <<std::endl;
    }
}       

const NxsSimpleNode * downpassSettingNodesAboveFields(const NxsSimpleNode * rootLeft, const std::vector<const NxsSimpleNode *> & preorderTraversal,std::vector<NdBlob *> &gAllocedBlobs) {
    if (gDebugging) {
        std::cerr << "Starting downpass" <<std::endl;
    }
    const NxsSimpleNode * centroidChild = 0L;
    int minCentroidScore = gNumLeaves;
    NdBlob * rlb = new NdBlob(true);
    rlb->SetIsParentsLeftChild(true);
    rootLeft->scratch = (void *) rlb;
    rlb->numLeavesAboveEdge = 1;
    rlb->SetNumActiveLeavesAbove(1);
    EdgeDecompInfo * rlEdgeInfo = rlb->GetActiveEdgeInfoPtr();
    rlEdgeInfo->GetClosestLeavesAboveRef().push_back(LeafPathElement(rootLeft, THIS_NODE_DIR));
    rlEdgeInfo->SetAboveInitialized(true);
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
            edgeInfo->GetClosestLeavesAboveRef().push_back(LeafPathElement(nd, THIS_NODE_DIR));
            edgeInfo->SetAboveInitialized(true);
        }
        else {
            nb->numLeavesAboveEdge = LeftBlob(*nd)->numLeavesAboveEdge + RightBlob(*nd)->numLeavesAboveEdge;
            nb->SetNumActiveLeavesAbove(nb->numLeavesAboveEdge);
            mergeShortLeavesLists(*nd, LEFT_OR_RIGHT, *LeftChild(*nd), *RightChild(*nd));
            int currSC = CentroidScore(*nd, gNumLeaves);
            if (currSC < minCentroidScore) {
                centroidChild = nd;
                minCentroidScore = currSC;
                if (gDebugging) {
                    std::cerr << "New Min Centroid Score = " << minCentroidScore << std::endl;
                }
            }

        }
    }
    NdBlob * rb = new NdBlob(false);
    rb->SetIsParentsLeftChild(false);
    gRoot->scratch = (void *) rb;
    rb->numLeavesAboveEdge = LeftBlob(*gRoot)->numLeavesAboveEdge + RightBlob(*gRoot)->numLeavesAboveEdge;
    rb->SetNumActiveLeavesAbove(rb->numLeavesAboveEdge);
    rb->SetActiveLeafDir(LEFT_OR_RIGHT);
    
    
    gRoot->SetTaxonIndex(gNumLeaves);
    mergeShortLeavesLists(*gRoot, LEFT_OR_RIGHT, *LeftChild(*gRoot), *RightChild(*gRoot));
    assert(rb->GetActiveEdgeInfoPtr()->GetCloseLeavesAbove()[0].GetLeaf() == rootLeft);
    rb->GetActiveEdgeInfoPtr()->GetClosestLeavesAboveRef()[0].SetScore(0); // modify the score so that the root-spanning path will be treated as a "broken" edge
    rb->GetActiveEdgeInfoPtr()->SetBelowInitialized(true); // the root's "below" field is also set (empty).
    if (gDebugging) {
        std::cerr << "Ending downpass" <<std::endl;
    }

    
    return centroidChild;
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
    try {
        const NxsSimpleNode * centroidChild = 0L;
        if (gNumLeaves <= gMaxSubProblemSize) {
            if (gOutputStream) {
                *gOutputStream << "Tree tree" << gCurrTreeIndex << " = " << ftd.GetNewick() << ";" << std::endl;
            }
        }
        else {
            if (gEmitIntermediates) {
                *gOutputStream << "Tree " << gCurrTreeIndex << " = [&U] ";
                writeNewickOfActive(*gOutputStream, Parent(*rootLeft));
                *gOutputStream << ";\n";
            }
            centroidChild = downpassSettingNodesAboveFields(rootLeft, preorderTraversal, gAllocedBlobs);
            uppassSettingNodesBelowFields(preorderTraversal);
            assert(centroidChild != 0);
            NxsString namePrefix;
            namePrefix += gCurrTreeIndex;
    
            if (gDebugging) {
                debugPrintTree(std::cerr, gRoot);
            }
            debugCheckActiveFlags(gRoot, gNumLeaves);
            decomposeAroundCentroidChild(centroidChild, centroidChild, gNumLeaves, namePrefix, gRoot);
        }

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
    out << "    -x output intermediate trees.\n\n";
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
        else if (filepath[1] == 'x')
            gEmitIntermediates = true;
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

    int rc = readInput(*inpStream, f, filename);

    if (rc == 0 && gOutputStream != 0L) {
        *gOutputStream << "END;\n";
    }
    return rc;
}


