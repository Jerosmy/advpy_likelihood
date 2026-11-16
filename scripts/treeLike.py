#! /opt/anaconda3/bin/python

class nodes:
    def __init__(self, seq:str = str(), parent:str = str(), child:list = [], edgeIn:dict = {}, edgeOut:dict = {}):
        self.seq = seq       #character list, store the sequence of the object
        self.parent = parent   #char/int (to be decided), store the ID of the parent of the object
        self.child = child     #list of children ids, store a list of children 
        self.edgeIn = edgeIn    #dictionary {parent1: length1}, Gives edges coming to this node
        self.edgeOut = edgeOut   #dictionary {child1: length1, child2: length2]}, lists edges coming out of this node