#!/opt/anaconda3/bin/python

import numpy as np
import scipy
import math
from benchmark.timelog import timelog

class Tree: 
    def __init__(self):
        self.descendant={} #clé: parent | terme: enfants (un parent plusieurs enfant) - attribut 1
        self.ascendant={} #clé: enfant | terme: noeud parent (un descendant, un parent) - attribut 2
        self.branches={} #clé: enfant | terme: longueur de la branche - attribut 3
        self.sequences={} #clé: feuille | terme: séquence associée - attribut 4
        self.root=None #Noeuds sans parent = racien - attribut 5
        self.seqlength=None #longueur des séquences
        self.aligned=False  #vériier que toutes les séquences sont de mêmes longueurs

    @timelog
    def __str__(self):
        lines=[]
        for nom, seq in self.sequences.items():
            lines.append(f"> {nom}: {seq}")
        return "\n".join(lines)

    @timelog
    def __len__(self):
        if self.seqlength is None: 
            raise ValueError("sequences have not been loaded, or sequences of different lengths")
        return self.seqlength

    #Méthode 1: charger les données de parentés entre noeuds.
    @timelog
    def load_parente(self, nodes_path):
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
    @timelog
    def load_branches(self, branches_path):
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
    @timelog
    def load_sequences(self,sequences_path):
        with open(sequences_path,"r") as seq:
            content=seq.read().strip().split("\n")
            for ns in content:
                name, sequence=ns.split()
                self.sequences[name]=sequence
        self.check_lengths()

    
    #Méthode 4: définir l'ordre des calcules pour éviter de combiner des noeuds non résolus
    @timelog
    def define_order(self):
        order=[self.root]
        for nod in order:
            if nod in self.descendant:
                order.extend(self.descendant[nod])
        order.reverse()
        return order 

    @timelog
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

    @timelog
    def total_branch_lenght(self):
        return sum([bl for bl in self.branches.values()])



class Likelihood:
    @timelog
    def __init__(self,mu,tree,):
        self.mu=None
        self.Q=None
        self.transition_matrix(mu)
        self.onehotvectors={"A":[1,0,0,0],"C":[0,1,0,0],"G":[0,0,1,0],"T":[0,0,0,1]}
        self.pi=np.array([1/4,1/4,1/4,1/4])
        self.tree=tree
    
    #Méthode 1: construire la matrice de transition Q à partir de mu
    @timelog
    def transition_matrix(self, mu):
        self.mu=mu
        self.Q=np.array([[-3*mu,mu,mu,mu],[mu,-3*mu,mu,mu],[mu,mu,-3*mu,mu],[mu,mu,mu,-3*mu]])

    @timelog
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

    @timelog # Bottlneck (takes 10^3 more time to run)
    def prob_sequence(self):
        nbsites=len(self.tree)
        phyloprob=[]
        for site in range(nbsites):
            phyloprob.append(self.prob_site(site))
        return sum(phyloprob)
        #on pourrait ajouter la liste des résultats par site. 
        
    @timelog
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