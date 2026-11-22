#!/opt/anaconda3/bin/python

import numpy as np
import scipy
import math
from benchmark.timelog import timelog

class Tree: 
    def __init__(self):
        """ 
        Creates a Tree class object
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
        self.order = None # Save order of branches instead of always recomputing

    #@timelog
    def __str__(self):
        lines=[]
        for nom, seq in self.sequences.items():
            lines.append(f"> {nom}: {seq}")
        return "\n".join(lines)

    #@timelog
    def __len__(self):
        if self.seqlength is None: 
            raise ValueError("sequences have not been loaded, or sequences of different lengths")
        return self.seqlength

    #Méthode 1: charger les données de parentés entre noeuds.
    #@timelog
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
                
    #Méthode 2: charger les données correspondant à la longueur entre les noeuds. 
    #@timelog
    def load_branches(self, branches_path):
        """ 
        Loads branch lengths in the Tree object
        Args:
            branches_path (str): Path to the branch length .dat file.Values should be in the same order as the node file and be separated by a ","
        """
        #Cette manière de faire correspond avec l'exemple du cours sur le pdf. 
        #Parmi les enfants, si pas retrouvé chez les parents alors feuille terminale
        leaf=[nods for nods in self.ascendant if nods not in self.descendant]
        
        #Parmi les parents si retrouvés chez les enfants alors noeuds interns
        internalnodes=[nods for nods in self.descendant if nods in self.ascendant]

        #Liste réodronnant les enfants de sorte à être alignés avec les position de longueurs des branches. 
        ordnods=leaf+internalnodes
        
        with open(branches_path, "r") as brlen:
            longueurs=[float(dis) for dis in brlen.read().replace("\n","").split(",")]
            #for x in brlen.read().replace("\n","").split(","):
                #longeurs.append(float(x))
            for bl in range(len(longueurs)):
               self.branches[ordnods[bl]]=longueurs[bl]

    #Méthode 3: charger les séquences associées aux feuilles de l'arbre. 
    #@timelog
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
    #@timelog
    def define_order(self):
        if self.order is None:
            order = [self.root]
            for nod in order:
                if nod in self.descendant:
                    order.extend(self.descendant[nod])
            order.reverse()
            self.order = order
        
        return self.order

    #@timelog
    #Méthode 5: vérification des longueurs
    def check_lengths(self):
        """ 
        verifies alignement of sequences 
        """
        length=[len(seq) for seq in self.sequences.values()]
        if len(set(length))==1:
            self.aligned=True
            self.seqlength=length[0]
        else:
            self.aligned=False
            self.seqlength=None
            print(f"Warning: different sequences lengths : {length}")

    #@timelog
    def total_branch_lenght(self):
        return sum([bl for bl in self.branches.values()])



class Likelihood:
    #@timelog
    def __init__(self,mu,tree,):
        """ 
        Creates a Likelihood object
        """
        self.mu=None
        self.Q=None
        self.transition_matrix(mu)
        self.onehotvectors={"A":[1,0,0,0],"C":[0,1,0,0],"G":[0,0,1,0],"T":[0,0,0,1]}
        self.pi=np.array([1/4,1/4,1/4,1/4])
        self.tree=tree
        self.prob_matrix = {}
    
    #Méthode 1: construire la matrice de transition Q à partir de mu
    ##@timelog
    def transition_matrix(self, mu):
        """ 
        Constructs and attributes the transition matrix
        Args:
            mu (float): a transition rate
        Reuturns:
            A np array stored in the Likelihood.mu attribute
        """
        self.mu=mu
        self.Q=np.array([[-3*mu,mu,mu,mu],[mu,-3*mu,mu,mu],[mu,mu,-3*mu,mu],[mu,mu,mu,-3*mu]])

    def prob_mat(self):
        """ 
        Precomputes branch-wise probability matrixes
        Args :
            method on a Tree object
        Returns:
            A precomputed transition matrix
        """
        for a in self.tree.define_order():
            if a != self.tree.root:
                t=self.tree.branches[a]
                proba_matrix=scipy.linalg.expm(self.Q*t)
                self.prob_matrix[a] = proba_matrix
            else :
                pass

        

    
    ##@timelog
    def prob_site(self, site): 
        """ 
        
        """
        descendant=self.tree.descendant
        root=self.tree.root
        probs={}
        for a in self.tree.order:
            if a !=root:
                #t=self.tree.branches[a]
                #proba_matrix=scipy.linalg.expm(self.Q*t)
                proba_matrix = self.prob_matrix[a]
                if a not in descendant:
                    ohv=self.onehotvectors[self.tree.sequences[a][site]]
                    p=np.matmul(proba_matrix,ohv)
                    probs[a]=p
                else:
                    k1,k2=descendant[a]
                    pk1=probs[k1]
                    pk2=probs[k2]
                    vec = np.asarray(pk1) * np.asarray(pk2)
                    p=np.matmul(proba_matrix,vec)
                    probs[a]=p
            else: 
                k1,k2=descendant[a]
                pk1=probs[k1]
                pk2=probs[k2]
                vec = np.asarray(pk1) * np.asarray(pk2)
                probs[a]=vec
        npmm=np.matmul(probs[root],self.pi)
        return math.log(npmm)

    #@timelog # Bottlneck (takes 10^3 more time to run)
    def prob_sequence(self):
        
        return sum([self.prob_site(site) for site in range(len(self.tree))])
        #nbsites=len(self.tree)
        #phyloprob=[]
        #for site in range(nbsites):
        #    phyloprob.append(self.prob_site(site))
        #return sum(phyloprob)
        #on pourrait ajouter la liste des résultats par site. 
        
    ##@timelog
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