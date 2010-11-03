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

enum TreeSweepDirection {
    LEFT_OR_RIGHT = 0,
    LEFT_OR_BELOW = 1,
    RIGHT_OR_BELOW = 2
};
enum TreeDirection {
    LEFT_DIR = 0,
    RIGHT_DIR = 1,
    PARENT_DIR = 2,
    PAR_LEFT_LEFT = 3,  // for dealing with the root branch being broken
    PAR_LEFT_RIGHT = 4, // for dealing with the root branch being broken
    PAR_RIGHT_LEFT = 5, // for dealing with the root branch being broken
    PAR_RIGHT_RIGHT = 6, // for dealing with the root branch being broken
    THIS_NODE_DIR = 7
};

class LeafPathElement {
    public:
    LeafPathElement(const NxsSimpleNode *nd, TreeDirection d) {
        assert(nd);
        assert(nd->IsTip());
        leaf = nd;
        iScore = 1;
        dirToNext = d;
    }
    
    bool operator<(const LeafPathElement & other) const {
        return this->iScore < other.iScore;
    }
    
    private:
        const NxsSimpleNode * leaf;
        long iScore;
        TreeDirection dirToNext;
};

class NdDecompInfo {
    public:
        std::list<LeafPathElement> closestLeavesAbove;
        std::list<LeafPathElement> closestLeavesLeftOrBelow;
        std::list<LeafPathElement> closestLeavesRightOrBelow;
};

class NdBlob {
    public:
    
        NdBlob() {
            Reset();
        }
        
        long GetNumLeavesBelow() const {
            return gNumNodes - this->numLeavesAboveEdge;
        }
        void Reset() {
            this->numLeavesAboveEdge = -1;
            assert(ndInfoStack.empty());
            ndInfoStack.push(new NdDecompInfo());
            currNdInfo = ndInfoStack.top();
            //this->numLeavesBelowEdge = -1;
        }
        
        long numLeavesAboveEdge;
        std::stack<NdDecompInfo*> ndInfoStack;
        NdDecompInfo * currNdInfo;
        
//        long numLeavesBelowEdge;
};
std::vector<NdBlob> gAllNodeBlobs;
unsigned gMaxNumLeafPaths = 0;

// user-controlled options as global variables...
bool gVerbose = false;
long gLeafSetIntersectionSize = 50;
long gMaxSubProblemSize = 250;
bool gUseEdgeLengths = false;
bool gIntercalate = false;

std::ostream * gOutputStream = &(std::cout);


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



void mergePathElementLists(std::list<LeafPathElement> &peList,
                           TreeDirection firDir, 
                           const std::list<LeafPathElement> & firSource,
                           TreeDirection secDir, 
                           const std::list<LeafPathElement> & secSource
                           ) {
}

inline void mergeShortLeavesLists(const NxsSimpleNode &nd, TreeSweepDirection dir, const NxsSimpleNode &o, const NxsSimpleNode &o2) {
    if (gIntercalate) {
        assert(false);
        throw 5;
    }
    else {
        NdDecompInfo * ndInfo = ((NdBlob*)nd.scratch)->currNdInfo;
        assert(ndInfo);
        if (dir == LEFT_OR_RIGHT) {
            std::list<LeafPathElement> &peList = ndInfo->closestLeavesAbove;
            assert(peList.empty());
            const NxsSimpleNode * leftNd = &o;
            const NxsSimpleNode * rightNd = &o2;
            mergePathElementLists(peList, LEFT_DIR, ((NdBlob*)leftNd->scratch)->currNdInfo->closestLeavesAbove,
                                          RIGHT_DIR, ((NdBlob*)rightNd->scratch)->currNdInfo->closestLeavesAbove);
        }
        else {
            const NxsSimpleNode * belowNd = &o2;
            NdDecompInfo * parInfo = ((NdBlob*)belowNd->scratch)->currNdInfo;
            NdDecompInfo * oInfo = ((NdBlob*)o.scratch)->currNdInfo;
            bool currIsLeft = (&nd == LeftChild(*belowNd));
            const std::list<LeafPathElement> & parSource = (currIsLeft ? parInfo->closestLeavesRightOrBelow : parInfo->closestLeavesLeftOrBelow);
            if (dir == LEFT_OR_BELOW) {
                std::list<LeafPathElement> &peList = ndInfo->closestLeavesLeftOrBelow;
                const NxsSimpleNode * leftNd = &o;
                mergePathElementLists(peList, LEFT_DIR, oInfo->closestLeavesAbove,
                                              PARENT_DIR, parSource);
                
            }
            else {
                assert(dir == RIGHT_OR_BELOW);
                std::list<LeafPathElement> &peList = ndInfo->closestLeavesRightOrBelow;
                const NxsSimpleNode * rightNd = &o;
                mergePathElementLists(peList, RIGHT_DIR, oInfo->closestLeavesAbove,
                                              PARENT_DIR, parSource);
            }
        }
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
    


    NxsSimpleTree t(ftd, 0, 0.0);
    NxsSimpleNode * iniRoot = const_cast<NxsSimpleNode *>(t.GetRootConst());
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
        bogusNode.AddChild(iniRoot);
        gRoot = & bogusNode;
    }
    
    std::vector<const NxsSimpleNode *> preorderTraversal;
    gRoot->AddSelfAndDesToPreorder(preorderTraversal);

    gNumNodes = (long) preorderTraversal.size();
    if ((gNumNodes % 2) == 0) {
        gErrorMessage << "Tree " << gCurrTreeIndex + 1 << " has an even number of nodes!";
        throw NxsException(gErrorMessage);
    }
    gNumLeaves = (long)((gNumNodes + 1 ) / 2);
    if (gNumLeaves < gLeafSetIntersectionSize) {
        gErrorMessage << "LeafSetIntersectionSize = " << gLeafSetIntersectionSize << ", but there are only " << gNumLeaves << " leaves in tree " << gCurrTreeIndex + 1  << '.';
        throw NxsException(gErrorMessage);
    }
    
    
    std::vector<NxsSimpleNode *> bogusChildren = bogusNode.GetChildren();
    try {
        
        if (gNumLeaves <= gMaxSubProblemSize) {
            if (gOutputStream) {
                *gOutputStream << "Tree tree" << gCurrTreeIndex << " = " << ftd.GetNewick() << ";" << std::endl;
            }
        } 
        else {
            gAllNodeBlobs.resize((unsigned) gNumNodes);
            std::vector<NdBlob>::iterator nbIt = gAllNodeBlobs.begin();
            
            std::vector<const NxsSimpleNode *>::const_reverse_iterator postIt = preorderTraversal.rbegin();
            
            for (; postIt != preorderTraversal.rend(); ++postIt) {
                const NxsSimpleNode * nd = *postIt;
                NdBlob * nb = &(*nbIt++);
                nd->scratch = (void *) nb;
                nb->Reset();
                if (nd->IsTip()) {
                    nb->numLeavesAboveEdge = 1;
                    NdDecompInfo * ndInfo = nb->currNdInfo;
                    ndInfo->closestLeavesAbove.push_back(LeafPathElement(nd, THIS_NODE_DIR));
                    
                }
                else {
                    nb->numLeavesAboveEdge = LeftBlob(*nd)->numLeavesAboveEdge + RightBlob(*nd)->numLeavesAboveEdge;
                    mergeShortLeavesLists(*nd, LEFT_OR_RIGHT, *LeftChild(*nd), *RightChild(*nd));
                    
                }
                //std::cerr << "nd = " << (long) nd << " blob->numLeavesAboveEdge = " << nb->numLeavesAboveEdge << '\n';
            }
            
            if (((NdBlob *)gRoot->scratch)->numLeavesAboveEdge != gNumLeaves) {
                gErrorMessage << "Expected " << gNumLeaves << " leaves at the root, but only found " << ((NdBlob *)gRoot->scratch)->numLeavesAboveEdge;
                throw NxsException(gErrorMessage);
            }
            
            std::vector<const NxsSimpleNode *>::const_iterator preIt = preorderTraversal.begin();
            
            
            
            
        }        
        
        
    }
    catch (...) {
         for (std::vector<NxsSimpleNode *>::iterator cIt = bogusChildren.begin(); cIt != bogusChildren.end(); ++cIt)
            bogusNode.RemoveChild(*cIt);
        throw;
        }

    for (std::vector<NxsSimpleNode *>::iterator cIt = bogusChildren.begin(); cIt != bogusChildren.end(); ++cIt)
        bogusNode.RemoveChild(*cIt);
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
	out << "    -h on the command line shows this help message.\n\n";
	out << "    -v verbose output to stdandard error.\n\n";
	out << "    -e use EDGE_LENGTHS in the shortest leaf calculations.\n\n";
	out << "    -i<non-negative integer> specifies MIN_LEAF_INTERSECTION_SIZE\n";
	out << "    -m<non-negative integer> specifies MAX_SUBTREE_SIZE\n";
	out << "    -f<format> specifies the input file format expected:\n";
	out << "            -fnexus     NEXUS (this is also the default)\n";
	out << "            -fphyliptree    Phylip tree file (simple newick string)\n";
	out << "            -frelaxedphyliptree Relaxed Phylip name restrictions\n";
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
    int rc =  readInput(*inpStream, f, filename);
    if (rc == 0 && gOutputStream != 0L) {
        *gOutputStream << "END;\n";
    }
    return rc;
}


