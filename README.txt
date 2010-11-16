centroid_decomp  by Mark T. Holder

################################################################################
Building
################################################################################
1. Download the source code for NCL version 2.1  I recommend using svn to 
    get the latest version of the code:
 
   $ svn checkout https://ncl.svn.sourceforge.net/svnroot/ncl/branches/v2.1 ncl2.1

2. Build and install NCL.

   $ sh bootstrap.sh
   $ ./configure --prefix=/my/favorite/install/dir
   $ make 
   $ make check
   $ make install
   $ make installcheck

3. Set the environmental variable NCL_INSTALL_DIR to the directory used at the
    prefix to NCL's configure script (if you did not use a --prefix argument, 
    then the default is /usr/local):
 
    $ export NCL_INSTALL_DIR=/my/favorite/install/dir

4. Use make in the src directory:

    $ cd src
    $ make
    $ make check

################################################################################
Usage
################################################################################
Reads a file with a tree description and produces a decomposition of the tree
into smaller trees.  Subject to the following constraints:
    1. No tree will have more than MAX_SUBTREE_SIZE leaves
    2. Trees will be written in hierarchichal groups (with their names
           indicating their group membership)
    3. Within each group, each tree will contain a common subset of leaves.
           The size of this subset will be at least MIN_LEAF_INTERSECTION_SIZE.
    4. All leaves in the original tree will be present in at least one output tree.

The input tree must be fully-resolved (binary).

Decomposition is produced by recursively breaking the tree around centroid
    edge and striving to include a group of MIN_LEAF_INTERSECTION_SIZE/4 leaves
    from each induced subtree into the smaller tree (if some subtrees are 
    smaller than MIN_LEAF_INTERSECTION_SIZE/4 then more leaves will be included
    from their "sister" subtree).

At each step, leaves closest to the decomposition edge are retained.  This 
    distance can be topological or based on EDGE_LENGTHS (-e option below; in 
    this case all edges in the tree must have edge lengths).

Command-line flags:

    -h on the command line shows this help message.

    -v verbose output to stdandard error.

    -e use EDGE_LENGTHS in the shortest leaf calculations.

    -i<non-negative integer> specifies MIN_LEAF_INTERSECTION_SIZE
    -m<non-negative integer> specifies MAX_SUBTREE_SIZE
    -f<format> specifies the input file format expected:
            -fnexus     NEXUS (this is also the default)
            -fphylip    Phylip tree file (simple newick string)
            -frelaxedphylip Relaxed Phylip name restrictions


Example:
    ./centroid_decomp -frelaxedphyliptree myfile.phy  
    

################################################################################
Running times 
################################################################################
Timed on a Mac laptop with:
    Intel Core i7 - 2.66 GHz processor
    8G RAM
    ssd
    
Some random tree shapes:
   27,643 taxa ->  1.248 seconds
   55,286 taxa ->  2.530 seconds
  110,572 taxa ->  5.128 seconds
     
Balanced trees:
     512 taxa ->  0.021 seconds
   1,024 taxa ->  0.045 seconds
   2,048 taxa ->  0.087 seconds
   4,096 taxa ->  0.173 seconds
   8,192 taxa ->  0.349 seconds
  16,384 taxa ->  0.695 seconds
  32,768 taxa ->  1.409 seconds
  65,536 taxa ->  2.866 seconds
 131,072 taxa ->  5.853 seconds
 262,144 taxa -> 12.335 seconds
 524,288 taxa -> 27.057 seconds

Pectinate trees:
     512 taxa ->  0.024 seconds
   1,024 taxa ->  0.051 seconds
   2,048 taxa ->  0.098 seconds
   4,096 taxa ->  0.203 seconds
   8,192 taxa ->  0.396 seconds
  16,384 taxa ->  0.817 seconds
  32,768 taxa ->  1.678 seconds
  65,536 taxa ->  3.450 seconds
 131,072 taxa ->  7.173 seconds
 262,144 taxa -> 15.473 seconds

################################################################################
License
################################################################################
See Licensce.txt for the text of the BSD license used for this software.
