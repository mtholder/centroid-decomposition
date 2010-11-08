#if ! defined(CENTROID_DECOMP_HPP)
#define CENTROID_DECOMP_HPP

#include "ncl/nxsmultiformat.h"


extern long gNumNodes;
extern long gNumLeaves;



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

        NdBlob(bool parentsLeftC)
            :numActiveLeavesAboveEdge(-1),
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
        void SetNumActiveLeavesAbove(long n){
            this->numActiveAboveStack.push(this->numActiveLeavesAboveEdge);
            this->numActiveLeavesAboveEdge = n;
        }

        void SetActiveEdgeInfoPtr(EdgeDecompInfo *n) {
            this->edgeInfoStack.push(this->activeEdgeInfo);
            this->activeEdgeInfo = n;
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
        

        void SetActiveLeafDir(TreeSweepDirection n) {
            this->activeLeafDirStack.push(this->activeLeafDir);
            this->activeLeafDir = n;
        }
        
        TreeSweepDirection GetActiveLeafDir() const {
            return this->activeLeafDir;
        }

        TreeSweepDirection GetActiveAndFocalLeafDir() const {
            int faldi = ((int) this->activeLeafDir) | ((int) this->focalEdgeDir) ;
            return TreeSweepDirection(faldi);
        }

        void SetFocalEdgeDir(TreeDirection n) {
            this->focalEdgeDirStack.push(this->focalEdgeDir);
            this->focalEdgeDir = n;
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
        
        
        long numLeavesAboveEdge;
        int ndDirWRTParent;
        LPECollection lpeScratch1;
        LPECollection lpeScratch2;


    private:

        EdgeDecompInfo fullEdgeInfo;

        long numActiveLeavesAboveEdge;
        TreeSweepDirection activeLeafDir;
        EdgeDecompInfo * activeEdgeInfo;
        TreeDirection focalEdgeDir; // direction to move to find the focal edge (e.g. the centroid edge) This will be BELOW_DIR for the node is attached to the edge

        std::stack<long> numActiveAboveStack;
        std::stack<EdgeDecompInfo*> edgeInfoStack;
        std::stack<TreeSweepDirection> activeLeafDirStack;
        std::stack<TreeDirection> focalEdgeDirStack;

        bool isParentsLeftChild;

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



inline void writeLeafSet(std::ostream &o, const LPECollection &ls) {
    o << "{";
    for (LPECollectionConstIt i = ls.begin(); i != ls.end(); ++i) {
        if (i != ls.begin())
            o << ", ";
        const NxsSimpleNode &nd = *(i->leaf);
        o << nd.GetTaxonIndex();
    }
    o << "}";
}



#endif  /*! defined(CENTROID_DECOMP_HPP) */

