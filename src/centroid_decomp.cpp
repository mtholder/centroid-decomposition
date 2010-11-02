#include "ncl/nxsmultiformat.h"
#include <iostream>
#include <fstream>
#include <string>


void printHelp(std::ostream & out);


bool gVerbose = false;
std::string gVersionString("centroid_decomp version 0.0.1");
std::string gErrMsgPrefix("centroid_decomp: ");



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
	out << "            -fphylip    Phylip tree file (simple newick string)\n";
	out << "            -frelaxedphylip Relaxed Phylip name restrictions\n";
}


int main(int argc, char * argv[]) {
	NxsReader::setNCLCatchesSignals(true);
	MultiFormatReader::DataFormatType f(MultiFormatReader::NEXUS_FORMAT);
    std::string formatName("NEXUS");
    long leafSetIntersectionSize = 50;
    long maxSubProblemSize = 250;
    bool useEdgeLengths = false;
    
	for (int i = 1; i < argc; ++i) {
		const char * filepath = argv[i];
		const unsigned slen = strlen(filepath);
		if (slen < 2 || filepath[0] != '-')
			continue;
		if (filepath[1] == 'h') {
			printHelp(std::cout);
			return 1;
		}
		else if (filepath[1] == 'v')
			gVerbose = true;
		else if (filepath[1] == 'e')
			useEdgeLengths = true;
        else if (filepath[1] == 'i') {
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &leafSetIntersectionSize))) {
				std::cerr << gErrMsgPrefix << "Expecting an integer -i\n" << std::endl;
				return 3;
			}
		}
        else if (filepath[1] == 'm') {
			if ((slen == 2) || (!NxsString::to_long(filepath + 2, &maxSubProblemSize))) {
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

    if (leafSetIntersectionSize < 4) {
        std::cerr << gErrMsgPrefix << "Minimum size of the leaf intersection set must 4" << std::endl;
        return 6;
    }
    if (maxSubProblemSize < leafSetIntersectionSize + 3) {
        std::cerr << "Max Sub problem size must be at least 3 larger that the minimum size of the leaf intersection set" << std::endl;
        return 7;
    }


	bool filefound = false;
	std::string filename;
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
	    if (useEdgeLengths) {
	        std::cerr << " -e";
	    }
	    std::cerr << " -i" << leafSetIntersectionSize;
	    std::cerr << " -m" << maxSubProblemSize;
	    std::cerr << " -f" << formatName;
	    if (filefound) {
	        std::cerr << " \"" << filename  << "\"";
	    }
	    std::cerr << std::endl;
	}
	if (gVerbose && !filefound) {
	    std::cerr << gErrMsgPrefix << "Reading from stdin..." << std::endl;
	}

	
	
	return 0;
}


