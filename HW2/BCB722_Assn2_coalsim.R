# BCB 722 Topics in Population Genetics
# Assignment #2: Coalescent Simulations


# Preparing to Simulate the Coalescent:
# As discussed in the assignment document, we will need two specialized data structures.
# For the first, we will define a simple class to represent the nodes at the current
# level in our tree (i.e. one horizontal slice of time).
# We only need to keep track of two things for the purposes of this simulation:
#    1) for each node in the tree, the set of all individuals in our sample
#       that are descendants of the node
#    2) the time of this node (for our actual samples this is the present, or 0)
# Thus, each node will have a vector of descendants (i.e. samples, aka leaves in the
# tree, that are descendants of the node), and a time attribute as well
treeNode <- setClass("treeNode", slots = c(descendants = "numeric", time = "numeric"))

# similarly, we will keep track of the branches and their lengths; this information
# will make it very easy to add our mutations later. This class is identical to the
# treeNode class expect it will be used to keep track of branch lengths rather than
# node ages.
branch <- setClass("branch", slots = c(descendants = "numeric", branchLength = "numeric"))

# Make your own functions here ####
### Make the following three functions to be used within the coalescent simulator
pickCoalTime <- function(???){
  # Hint: What distribution describes the number of generations we have to go back in
  # time for our next coalescence event to occur 
  #coalGen <- ???
  return(coalGen)
}

### this should be a vector containing the two numbers: the indices of the two nodes that we would
# like to coalesce.
selectNodeIndicesToCoalesce <- function(???){  
  # n = number of nodes (gene copies) to coalesce
  #nodesToCoalesce <- ???
  return(nodesToCoalesce)
}

### returns the per-generation probability that there will be any coalescence event among our lineages 
calcCoalProb <- function(???, ???){
  # two input arguments are required
  #coalProb <- ???
  return(coalProb)
}

# Tree-generating function will call your functions above... ####
# Add your own code where you see triple hash marks "###"
generateTree <- function(n, N){
  # initialize our tree nodes to contain our leaves
  nodes <- list()
  for (i in 1:n){
    currNode <- treeNode(descendants=c(i), time=0)
    nodes[[i]] <- currNode
  }
  
  # initialize our branches (currently empty)
  branches <- list()
  
  # we want to keep track of how far back in time we have gone
  # we will use this variable to do that
  currentTimeElapsed <- 0
  
  lineageCount = n #keeping track of the number of lineages yet to coalesce
  while (lineageCount >= 2){
    ### calculate coalescence probability:
    coalProb <- calcCoalProb(???,???) ### Add arguments to call your own function here
    
    # here is one of the ingenious things about coalescent simulation:
    # we don't actually have to simulate each generation, and this is
    # because the coalescent tells us what the distribution of
    # coalescence times should be. all we have to do is draw from this
    # distribution to get a number of generations, then we skip back
    # to that time and continue our work. this makes the simulations
    # really fast! we do this below:
    #
    # pick a random number of generations we will have to wait for the
    # coalescence event to occur:    
    coalGen <- pickCoalTime(???) ### Your function (and argument(s)) in use here
    
    # pick a random pair of nodes from our tree that will coalesce
    indicesToCoalesce <- selectNodeIndicesToCoalesce(lineageCount)
    node1 <- nodes[[indicesToCoalesce[1]]]
    node2 <- nodes[[indicesToCoalesce[2]]]
    
    # we just jumped back coalGen generations, so we will add this to our
    # time counter
    currentTimeElapsed <- currentTimeElapsed + coalGen
    
    # now that node1 and node2 have coalesced in our tree, we remove them from
    # the set of nodes that have yet to coalesce. because this list is an ordered
    # data structure, we will remove the second of the two in our list. The first
    # will later be replaced by our new node
    nodes[[max(indicesToCoalesce)]] = NULL
    
    # now we must replace those nodes with their ancestor's node;
    # first we create this node, whose descendants are the union of those of
    # the two nodes that just coalesced
    # then we insert it into our list of nodes
    
    newNode =  ### Write your own code to create the newNode!
    
    nodes[[min(indicesToCoalesce)]] = newNode  # insert the newNode into our list of nodes
    # in essence we have added two new branches to the tree: one from node1 and
    # one from node2; these two branches connect to each other at currentTimeElapsed
    
    # eventually we will be throwing mutations down onto branches;
    # to do this we need to keep track of two things:
    #     1) the length of the branch (in number of generations)
    #     2) the set of nodes that are a descendant of this branch
    # the latter will tell us which individuals will harbor any mutation(s)
    # that will occur on this branch
    # first, let's get the length of the branch coming from node1:
    branchLength1 <- currentTimeElapsed - node1@time
    
    # let's record the branch leading from node1 to the node we just added
    newBranch1 <- branch(descendants=node1@descendants, branchLength=branchLength1)
    branches[[length(branches)+1]] <- newBranch1
    
    # do the same thing for the branch leading to node2
    branchLength2 <- currentTimeElapsed - node2@time
    newBranch2 <- branch(descendants=node2@descendants, branchLength=branchLength2)
    branches[[length(branches)+1]] <- newBranch2
    
    lineageCount = lineageCount - 1 #one step closer to the common ancestor!
  }
  return(list("nodes" = nodes, "branches" = branches))
}

### Create the  drawNumberOfMutationsOnBranch function. 
# This will be used in the 'throwDownMutations()' function below.
drawNumberOfMutationsOnBranch <- function(???, ???)
{
  # have to pick an appropriately distributed random number of mutations for
  # a given branch
  return(???)
}

# Mutations ####
# Notice that this function takes our branches list as input but doesn't
# use any information from our treeNode objects. The latter was used only
# to help us build the former. The branches list now has all of the
# information to throw down mutations.
throwDownMutations <- function(branches, n, mu, L)
{
    # initialize our matrix to n empty strings
    seqs = c()
    for (i in 1:n){
        seqs = c(seqs, "")
    }

    totalMuts <- 0
    # now we are ready to throw down our mutations. there are a few ways that we could
    # go about it, but the way I have written it below is to iterate through each
    # branch of the tree, pick an appropriate random number of mutations for that branch
    # (you will code that part), and then for each chromosome in our sample write a "1"
    # if it si a descendant of the branch where the mutation occurred and "0" otherwise.
    for (b in 1:length(branches)){
        currMuts = drawNumberOfMutationsOnBranch(branches[[b]]@branchLength, mu*L)
        mut <- 0
        while (mut < currMuts){
            for (i in 1:n){
                if (i %in% branches[[b]]@descendants){ #did individual i descend from branch b?
                    seqs[i] <- paste(seqs[i], "1", sep = "")
                }
                else{
                    seqs[i] <- paste(seqs[i], "0", sep = "")
                }
            }
            mut <- mut + 1
        }
        totalMuts <- totalMuts + currMuts
    }

    # we created an alignment of ones and zeros using string, with was nice and simple
    # but some of our downstream work might actually be easier to convert this into a
    # matrix of integers. so we do this below by initialize an appropriately sized
    # matrix of zeros and filling in the values from our sequences (converting to
    # integers as we go).
    seqMatrix <- matrix(0L, nrow=n, ncol=totalMuts)
    for (i in 1:n){
        j = 1
        while (j <= totalMuts){
            seqMatrix[i,j] <- as.integer(substr(seqs[i], j, j))
            j <- j + 1
        }
    }
    
    # in a non-recombining chromosome the order of our mutations should be
    # completely random, so we should shuffle the order of columns in our
    # mutation matrix (if we have any mutations, that is)
    if (totalMuts > 0){
        seqMatrix <- seqMatrix[,sample(1:totalMuts)]
    }
    return(seqMatrix)
}

### Create the following getTotalTreeHeight function: ####
# HINT: this can easily be obtained from one of our two tree data structure lists
getTotalTreeHeight <- function(???)
{
  return(???)
}

# Coalescent Simulation Parameters -- you can mess with these a little if you want
# but if you make L and N too huge things might get slow. Your results might also 
# start to act a little weird if you make n very large.
L = 100000  # Locus length
N = 1e4  # Population Size
n = 20  # number of pairs of genes that can coalesce
mu = 1.2e-8  # per site mutation rate
# Simulate
tree <- generateTree(n, N)
tmrca <- getTotalTreeHeight(???)  ### Using the function you made above...
print(tmrca)
seqMatrix <- throwDownMutations(tree$branches, n, mu, L)
print(seqMatrix)

### Your Analysis here:
# Site frequency Spectrum
# Come up with a way to generate the SFS from the simulated seqMatrix we generated above
# The hist() function is sufficient to make the plot
