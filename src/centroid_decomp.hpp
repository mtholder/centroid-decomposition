#if ! defined(CENTROID_DECOMP_HPP)
#define CENTROID_DECOMP_HPP

#include "ncl/nxsmultiformat.h"


extern long gNumNodes;
extern long gNumLeaves;
extern bool gDebugging;


enum TreeSweepDirection {
    NO_DIR_BIT = 0,
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
            iScore = 0;
            dirToNext = d;
            indexInNext = -1;
            indexInFull = 0;
        }
    
        LeafPathElement(const LeafPathElement &nextEl, TreeDirection d, int indexAtAdj, int indInFullList) {
            leaf = nextEl.leaf;
            dirToNext = d;
            indexInNext = indexAtAdj;
            indexInFull = indInFullList;
            assert(nextEl.iScore >= 0);
            iScore = nextEl.iScore + 1;
        }
    
        bool operator<(const LeafPathElement & other) const {
            assert(this->iScore >= 0);
            assert(other.iScore >= 0);
            return this->iScore < other.iScore;
        }
        
        long GetScore() const {
            return this->iScore;
        }
        void SetScore(long sc) {
            this->iScore = sc;
        }
        const NxsSimpleNode * GetLeaf() const {
            return this->leaf;
        }
    public: // should be private
        TreeDirection dirToNext;
        int indexInNext;
        unsigned indexInFull;
    private:
        const NxsSimpleNode * leaf;
        long iScore;
};
typedef std::vector<LeafPathElement> LPECollection;
typedef LPECollection::const_iterator LPECollectionConstIt;
typedef LPECollection::iterator LPECollectionIt;
std::ostream & operator<<(std::ostream &, const LeafPathElement &);

class EdgeDecompInfo {
    public:
        const LPECollection & GetCloseLeavesAbove() const {
#           if ! defined(NDEBUG)
                if (!aboveInitialized) {
                    std::cerr << "\nGetClosestLeavesAbove valid assertion failing for EDI allocated at line " << this->line << " of file " << this->fn << std::endl;
                    assert(aboveInitialized);
                }
#           endif 
            return closestLeavesAbove;
        }
        const LPECollection & GetCloseLeavesBelow() const {
#           if ! defined(NDEBUG)
                if (!belowInitialized) {
                    std::cerr << "\nGetClosestLeavesBelow valid assertion failing for EDI allocated at line " << this->line << " of file " << this->fn << std::endl;
                    assert(belowInitialized);
                }
#           endif ! defined(NDEBUG)
            return closestLeavesBelow;
        }

        LPECollection & GetClosestLeavesAboveRef() {
            return closestLeavesAbove;
        }
        LPECollection & GetClosestLeavesBelowRef() {
            return closestLeavesBelow;
        }
        void SetClosestLeavesAbove(const LPECollection & o) {
            SetAboveInitialized(true);
            closestLeavesAbove = o;
        }
        void SetClosestLeavesBelow(const LPECollection & o) {
            SetBelowInitialized(true);
            closestLeavesBelow = o;
        }

        void SetAboveInitialized(bool v=true) {
            this->aboveInitialized = v;
        }
        void SetBelowInitialized(bool v=true) {
            this->belowInitialized = v;
            
        }
        bool IsAboveInitialized() const {
            return this->aboveInitialized;
        }
        bool IsBelowInitialized() const  {
            return this->belowInitialized;
        }
        
        
#       if ! defined(NDEBUG)
            EdgeDecompInfo(const char *filename, int lineInit)
                :fn(filename),
                line(lineInit),
#       else
            EdgeDecompInfo(const char *, int )
                :
#       endif
                aboveInitialized(false),
                belowInitialized(false)
                {}
            
    private:
#       if ! defined(NDEBUG)
            std::string fn ; //debugging purposes only -- the line at which the instance was created.
            int line ; //debugging purposes only -- the line at which the instance was created.
#       endif
        LPECollection closestLeavesAbove;
        LPECollection closestLeavesBelow;
        bool aboveInitialized;
        bool belowInitialized;
        
};

class NdBlob;

class NdBlobSettingStruct {
    public:
        NdBlobSettingStruct(NdBlob *b, 
                        TreeDirection focal,
                        long numActiveAboveEdge,
                        TreeSweepDirection activeLeafSweepDir,
                        EdgeDecompInfo *edi,
                        std::set<NdBlob *> &altered)
            :blob(b),
            focalDir(focal),
            numActiveAbove(numActiveAboveEdge),
            activeLeafDir(activeLeafSweepDir),
            edgeDecompInfo(edi),
            alteredBlobs(altered)
            {
        }
        void SetObject();

    public:
        NdBlob * blob;
        TreeDirection focalDir;
        long numActiveAbove;
        TreeSweepDirection activeLeafDir;
        EdgeDecompInfo * edgeDecompInfo;
        std::set<NdBlob *> & alteredBlobs;
};

inline unsigned CheckFullIndexForLeaf(const NxsSimpleNode *nd, const std::map<const NxsSimpleNode *, int> & m) {
    std::map<const NxsSimpleNode *, int>::const_iterator toInd = m.find(nd);
    if (toInd == m.end()) {
        return UINT_MAX;
    }
    return toInd->second;
}

class NdBlob {
    public:

        NdBlob(bool parentsLeftC)
            :fullEdgeInfo(__FILE__, __LINE__),
            numActiveLeavesAboveEdge(-1),
            activeLeafDir(NO_DIR_BIT),
            focalEdgeDir(BELOW_DIR),
            isParentsLeftChild(parentsLeftC){
            Reset();
        }

        long GetNumLeavesBelow() const {
            return gNumNodes - this->numLeavesAboveEdge;
        }
        void Reset();
        void popBlob();

        long GetNumLeavesAbove() const {
            return this->numLeavesAboveEdge;
        }
        long GetNumActiveLeavesAbove() const {
            return this->numActiveLeavesAboveEdge;
        }

        
        EdgeDecompInfo * GetActiveEdgeInfoPtr() {
            return this->activeEdgeInfo;
        }
        const EdgeDecompInfo * GetActiveEdgeInfoConstPtr() const {
            return this->activeEdgeInfo;
        }
        
        EdgeDecompInfo * GetFullEdgeInfoPtr() {
            return &this->fullEdgeInfo;
        }
        const EdgeDecompInfo * GetFullEdgeInfoConstPtr() const {
            return &this->fullEdgeInfo;
        }
        

        TreeSweepDirection GetActiveLeafDir() const {
            return this->activeLeafDir;
        }

        TreeSweepDirection GetActiveAndFocalLeafDir() const {
            int faldi = ((int) this->activeLeafDir) | ((int) this->focalEdgeDir) ;
            return TreeSweepDirection(faldi);
        }

        
        TreeDirection GetFocalEdgeDir() const {
            return this->focalEdgeDir;
        }
        bool IsParentsLeftChild() const {
            return isParentsLeftChild;
        }
        void SetIsParentsLeftChild(bool v) {
            isParentsLeftChild = v;
        }
        
        const std::map<const NxsSimpleNode *, int> & GetLeafToIndexMapUp() const {
            return leafToIndexMapUp;
        }
        const std::map<const NxsSimpleNode *, int> & GetLeafToIndexMapDown() const {
            return leafToIndexMapDown;
        }
        
        long numLeavesAboveEdge;
        int ndDirWRTParent;
        LPECollection lpeScratch1;
        LPECollection lpeScratch2;
        mutable std::map<const NxsSimpleNode *, int> leafToIndexMapUp;
        mutable std::map<const NxsSimpleNode *, int> leafToIndexMapDown;

        
        unsigned FindFullIndexForLPEDown(const LeafPathElement & lpe, TreeDirection d, int indexInAdjacent) const {
            unsigned i = CheckFullDownIndexForLeaf(lpe.GetLeaf());
            if (i == UINT_MAX) {
                i = fullEdgeInfo.GetClosestLeavesBelowRef().size();
                fullEdgeInfo.GetClosestLeavesBelowRef().push_back(LeafPathElement(lpe, d, indexInAdjacent, (int) i));
                leafToIndexMapDown[lpe.GetLeaf()] = i;
            }
            return i;
        }
        unsigned FindFullIndexForLPEUp(const LeafPathElement & lpe, TreeDirection d, int indexInAdjacent) const {
            unsigned i = CheckFullUpIndexForLeaf(lpe.GetLeaf());
            if (i == UINT_MAX) {
                i = fullEdgeInfo.GetClosestLeavesAboveRef().size();
                fullEdgeInfo.GetClosestLeavesAboveRef().push_back(LeafPathElement(lpe, d, indexInAdjacent, (int) i));
                leafToIndexMapUp[lpe.GetLeaf()] = i;
            }
            return i;
        }

    private:
        unsigned CheckFullDownIndexForLeaf(const NxsSimpleNode *nd) const {
            return CheckFullIndexForLeaf(nd, leafToIndexMapDown);
        }
        unsigned CheckFullUpIndexForLeaf(const NxsSimpleNode *nd) const {
            return CheckFullIndexForLeaf(nd, leafToIndexMapUp);
        }
        void SetActiveEdgeInfoPtr(EdgeDecompInfo *n) {
            this->edgeInfoStack.push(this->activeEdgeInfo);
            this->activeEdgeInfo = n;
        }
        void SetActiveLeafDir(TreeSweepDirection n);
        
        void SetFocalEdgeDir(TreeDirection n) {
            this->focalEdgeDirStack.push(this->focalEdgeDir);
            this->focalEdgeDir = n;
        }
        void SetNumActiveLeavesAbove(long n){
            this->numActiveAboveStack.push(this->numActiveLeavesAboveEdge);
            this->numActiveLeavesAboveEdge = n;
        }

        mutable EdgeDecompInfo fullEdgeInfo;

        long numActiveLeavesAboveEdge;
        TreeSweepDirection activeLeafDir;
        EdgeDecompInfo * activeEdgeInfo;
        TreeDirection focalEdgeDir; // direction to move to find the focal edge (e.g. the centroid edge) This will be BELOW_DIR for the node is attached to the edge

        std::stack<long> numActiveAboveStack;
        std::stack<EdgeDecompInfo*> edgeInfoStack;
        std::stack<TreeSweepDirection> activeLeafDirStack;
        std::stack<TreeDirection> focalEdgeDirStack;

        bool isParentsLeftChild;
        
        
        friend const NxsSimpleNode *  downpassSettingNodesAboveFields(const NxsSimpleNode * rootLeft, const std::vector<const NxsSimpleNode *> & preorderTraversal,std::vector<NdBlob *> &gAllocedBlobs);
        friend void uppassSettingNodesBelowFields(const std::vector<const NxsSimpleNode *> & preorderTraversal);
        friend class NdBlobSettingStruct;
//        long numLeavesBelowEdge;
};





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

inline const NxsSimpleNode * SibNode(const NxsSimpleNode &nd) {
    const NxsSimpleNode *p = nd.GetParent();
    assert(p);
    if (((NdBlob *)nd.scratch)->IsParentsLeftChild())
        return RightChild(*p);
    return LeftChild(*p);
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

inline NdBlob * CurrentBlob(const NxsSimpleNode &nd) {
    return (NdBlob *)nd.scratch;
}


inline void writeLeafSet(std::ostream &o, const LPECollection &ls) {
    o << "ls@ " << (long)&ls <<" {";
    for (LPECollectionConstIt i = ls.begin(); i != ls.end(); ++i) {
        if (i != ls.begin())
            o << ", ";
        const NxsSimpleNode &nd = *(i->GetLeaf());
        o << nd.GetTaxonIndex();
    }
    o << "}";
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
    else if (d == NO_DIR_BIT)
        o << "NO_DIR_BIT";
    else {
        assert(false);
        o << "Unknown direction";
        }
    return o;
}

inline void writeLeafPathElement(std::ostream &o, const LeafPathElement &ls) {
    o << "lpe{leaf=" << ls.GetLeaf()->GetTaxonIndex();
    o << ", dirToNext=" << ls.dirToNext;
    o << ", indexInNext=" << ls.indexInNext;
    o << ", iScore=" << ls.GetScore() << '}';
}
inline std::ostream & operator<<(std::ostream & o, const LeafPathElement &lpe) {
    writeLeafPathElement(o, lpe);
    return o;
}




inline void writeLeafPathElementVector(std::ostream &o, const LPECollection &ls) {
    o << "{";
    for (LPECollectionConstIt i = ls.begin(); i != ls.end(); ++i) {
        if (i != ls.begin())
            o << ", ";
        writeLeafPathElement(o, *i);
    }
    o << "}";
}

inline void NdBlob::SetActiveLeafDir(TreeSweepDirection n) {
    this->activeLeafDirStack.push(this->activeLeafDir);
    this->activeLeafDir = n;
}


void decomposeAroundCentroidChild(const NxsSimpleNode *topCentroidChild,
                                  const NxsSimpleNode *bottomCentroidChild,
                                  long numActiveLeaves,
                                  NxsString namePrefix,
                                  const NxsSimpleNode * effectiveRoot);

#if ! defined(NDEBUG)
    void debugBlobPrint(std::ostream & o, const NdBlob * b);
#endif
void writeEdgeDecompInfo(std::ostream & o, const EdgeDecompInfo & edi);


inline void writeEdgeDecompInfo(std::ostream & o, const EdgeDecompInfo & edi) {
    o << "edi.CLA = ";
    if (edi.IsAboveInitialized()) 
        writeLeafSet(o, edi.GetCloseLeavesAbove());
    else
        o << "<uninit.>";
    o << "; edi.CLB = ";
    if (edi.IsBelowInitialized()) 
        writeLeafSet(o, edi.GetCloseLeavesBelow());
    else
        o << "<uninit.>";
};


#if ! defined(NDEBUG)
    inline void debugBlobPrint(std::ostream & o, const NdBlob * b) {
        o << "b.address = " << (long) b;
        o << " b#LvsAbove = " << b->GetNumLeavesAbove();
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
            o << " b.activeEdgeInfo@" << (long)b->GetActiveEdgeInfoConstPtr() << " = [";
            writeEdgeDecompInfo(o, *b->GetActiveEdgeInfoConstPtr());
            o << "]";
        }
    }
#endif
void mergePathElementLists(LPECollection &peList,
                           NdBlob * ndBlob,
                           TreeSweepDirection ,
                           TreeDirection firDir,
                           const LPECollection & firSource,
                           TreeDirection secDir,
                           const LPECollection & secSource
                           );

#endif  /*! defined(CENTROID_DECOMP_HPP) */

