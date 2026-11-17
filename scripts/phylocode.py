#!/opt/anaconda3/bin/python

import numpy as np
import scipy
import math


class Tree:

    def __init__(self):
        """ 
        Attributes:
            .descendant (dict): parent:children
            .ascendant (dict): child:parent node
            .branches (dict): child:branch length
            .sequences (dict): leave:sequence
            .root (str): root node, that has no parent, computed automatically
            .sequlength (int): length of sequences (should be aligned, i.e. all of saem len)
            .aligned (boolean): are all sequences of same length (i.e. aligned)?
        """ 
        self.descendant={} #clé: parent | terme: enfants (un parent plusieurs enfant) - attribut 1
        self.ascendant={} #clé: enfant | terme: noeud parent (un descendant, un parent) - attribut 2
        self.branches={} #clé: enfant | terme: longueur de la branche - attribut 3
        self.sequences={} #clé: feuille | terme: séquence associée - attribut 4
        self.root=None #Noeuds sans parent = racien - attribut 5
        self.seqlength=None #longueur des séquences
        self.aligned=False  #vériier que toutes les séquences sont de mêmes longueurs

    def __str__(self):
        lines=[]
        for nom, seq in self.sequences.items():
            lines.append(f"> {nom}: {seq}")
        return "\n".join(lines)

    def __len__(self):
        if self.seqlength is None: 
            raise ValueError("sequences have not been loaded, or sequences of different lengths")
        return self.seqlength

    #Method 1
    def load_parente(self, nodes_path):
        """ 
        Attirbutes parent-child relationships
        Args:
            nodes_path (str): system path to a per node child-parent relationship .dat table.

        Returns:
            Nothing, automatic completion of .ascendent and .descendent attributes of a Tree object
        """
        with open(nodes_path, "r") as nods:
            for line in nods:
                line = line.replace("\n","")
                if line =="":
                    continue
                ascendant, descendant=line.split(",")
                if ascendant not in self.descendant:
                    self.descendant[ascendant]=[]
                self.descendant[ascendant].append(descendant)
                self.ascendant[descendant] = ascendant
            for rac in self.descendant:
                if rac not in self.ascendant:
                    self.root=rac
                    break
                
    #Method 2
    def load_branches(self, branches_path:str = "../data/branchlength.dat"):
            """ 
            Loads branch lengths in the Tree object
            Args:
                branches_path (str): Path to the branch length .dat file.Values should be in the same order as the node file and be separated by a ","
            """
            with open(branches_path,"r") as brlen:
                branches_length=[float(bl) for bl in brlen.read().strip().split(",")]
                childs=list(self.ascendant)
                if len(branches_length)!=len(childs):
                    raise ValueError("Kinship not loaded or number of branches isn't equal to the number of child")
                else:
                    for CID in range(len(childs)): #CID: Child ID, position in the list
                        self.branches[childs[CID]]=branches_length[CID]

    #Method 3 
    def load_sequences(self,sequences_path):
        """ 
        Load the sequences of the organisms from a .msa file.

        Args:
            sequences_path (str): Path to the sequences in a msa format 
        """

        with open(sequences_path,"r") as seq:
            content=seq.read().strip().split("\n")
            for ns in content:
                name, sequence=ns.split()
                self.sequences[name]=sequence
        self.check_lengths()

    #Méthode 4: définir l'ordre des calcules pour éviter de combiner des noeuds non résolus
    def define_order(self):
        order=[self.root]
        for nod in order:
            if nod in self.descendant:
                order.extend(self.descendant[nod])
        order.reverse()
        return order 

    #Méthode 5: vérification des longueurs
    def check_lengths(self):
        length=[len(seq) for seq in self.sequences.values()]
        if len(set(length))==1:
            self.aligned=True
            self.seqlength=length[0]
        else:
            self.aligned=False
            self.seqlength=None
            print(f"Warning: different sequences lengths : {length}")


    def total_branch_lenght(self):
        return sum([bl for bl in self.branches.values()])



class Likelihood:
    def __init__(self,mu,tree,):
        self.mu=None
        self.Q=None
        self.transition_matrix(mu)
        self.onehotvectors={"A":[1,0,0,0],"C":[0,1,0,0],"G":[0,0,1,0],"T":[0,0,0,1]}
        self.pi=np.array([1/4,1/4,1/4,1/4])
        self.tree=tree
    #Méthode 1: construire la matrice de transition Q à partir de mu
    def transition_matrix(self, mu):
        self.mu=mu
        self.Q=np.array([[-3*mu,mu,mu,mu],[mu,-3*mu,mu,mu],[mu,mu,-3*mu,mu],[mu,mu,mu,-3*mu]])

    def prob_site(self, site): #revoir ce qu'on a envie de faire, tout calculer d'un coup
        descendant=self.tree.descendant
        root=self.tree.root
        probs={}
        for a in self.tree.define_order():
            if a !=root:
                t=self.tree.branches[a]
                if a not in descendant:
                    ohv=self.onehotvectors[self.tree.sequences[a][site]]
                    pm=scipy.linalg.expm(self.Q*t)
                    p=np.matmul(pm,ohv)
                    probs[a]=p
                else:
                    k1,k2=descendant[a]
                    pk1=probs[k1]
                    pk2=probs[k2]
                    vec=[float(pk1[i]*pk2[i]) for i in range(4)]
                    pm=scipy.linalg.expm(self.Q*t)
                    p=np.matmul(pm,vec)
                    probs[a]=p
            else: 
                k1,k2=descendant[a]
                pk1=probs[k1]
                pk2=probs[k2]
                vec=[float(pk1[i]*pk2[i]) for i in range(4)]
                probs[a]=vec
        npmm=np.matmul(probs[root],self.pi)
        return math.log(npmm)
    def prob_sequence(self):
        nbsites=len(self.tree)
        phyloprob=[]
        for site in range(nbsites):
            phyloprob.append(self.prob_site(site))
        return sum(phyloprob)
        #on pourrait ajouter la liste des résultats par site. 
        
    def variable_mu(self, mus):
        results=[]
        for mu in mus:
            self.transition_matrix(mu)
            probmuvar=self.prob_sequence()
            results.append(probmuvar)
            #on pourrait ajouter la liste des résultats par site. 
        return results



if __name__ == "__main__":
    Likelihood()
    Tree()