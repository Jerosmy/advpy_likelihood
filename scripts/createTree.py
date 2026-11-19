#!/opt/anaconda3/bin/python

from scripts import phylocode as pc
from benchmark.timelog import timelog

@timelog
def createTree(fpath="./data/table.dat", BRANCHPATH="./data/branchlength.dat", SEQPATH="./data/msa.dat"):
    #Object creation
    arbre=pc.Tree()

    #Loading the files
    arbre.load_parente(fpath)
    arbre.load_branches(BRANCHPATH)
    arbre.load_sequences(SEQPATH)

    #Define the order for branch length attribution
    arbre.define_order()

    ######################
    #Calculation part
    ######################

    ## New object creation 
    draft=pc.Likelihood(0.3, arbre)

    ##Computation
    draft_results= draft.prob_sequence()

    return draft_results


if __name__ == "__main__":
    print(createTree())