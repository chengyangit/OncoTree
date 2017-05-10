# -*- coding: utf-8 -*-

import sys
import os
import urllib2

import networkx as nx




class OncoTree():
    '''
    A OncoGxOne class to consume a tab delimited table of cancer type and construct tree structure class
    
    :Example:
    
    >>> from OncoTree import *
    >>> OT = OncoTree()
    >>> print dir(OT)  # Check all methods applicable to OncoTree class
    >>> ['OT', '__doc__', '__init__', '__module__', 'distance', 'homo', 'inclusion', 'url']
    >>> print OT.url
    >>> 'https://raw.githubusercontent.com/cBioPortal/oncotree/master/tumor_tree.txt'  # This is default cancer type input
    >>> print OT.abbrev['Leukemia']  # Check cancer type abbreviation. Note that the dictionary key has to be exact match. Since the cancer type terminology is mostly heterogeneous, please double check the cancer type metadata on url page specified by OT.url
    >>> 'LEUK'
    >>> print OT.OT
    >>> OT.OT # Return networkx clasee
    >>> networkx.classes.digraph.DiGraph object at 0x21f0450>
    >>> print OT.OT.successors('ALL')  # Check the successors of queried cancer type
    >>> ['TALL', 'BALL']  
    >>> print OT.OT.predecessors('ALL') # Check the predecessors of queried cancer type
    >>> ['LEUK']
    >>> print OT.inclusion('ALL')  # Check if input cancer type exist in OncoTree
    >>> True
    >>> print OT.homo('BALL', 'NSCLC')  # Check if two input cancer types are within same lineage
    >>> False
    >>> print OT.distance('BALL', 'NSCLC')  # Calculate the distance between two inpit cancer types 
    >>> 6 
    '''

 
    # fecth onco tree
    def __init__(self, file_path=False):
        self.url = "https://raw.githubusercontent.com/cBioPortal/oncotree/master/tumor_tree.txt"
        self.abbrev = {}
        
        # load the file.
        if not file_path:
    
            # fetch from inter-webs.
            req = urllib2.Request(self.url)
            response = urllib2.urlopen(req)
            the_page = response.read()
    
            # split into array.
            lines = the_page.strip().split("\n")
    
        else:
    
            # just open the file.
            lines = open(self.url, "rb")
    
        # create a graph.
        self.OT = nx.DiGraph()
    
        # add root node.
        self.OT.add_node("root", {"text": "root"})
        root = "root"
    
        # parse the file.
        line_cnt = 0
        for line in lines:
    
            # skip header.
            if line_cnt == 0:
                line_cnt += 1
                continue
    
            # tokenize.
            tokens = line.strip().split("\t")
            keys = tokens[0].split('(')[0].strip()
            values = tokens[0].split('(')[1][:-1]
            for col in range(5):
                if '(' in tokens[col]:
                    self.abbrev[tokens[col].split('(')[0].strip()] = tokens[col].split('(')[1][:-1]
            try:
                metamaintype = tokens[5]
            except IndexError:
                metamaintype = None
            try:
                metacolor = tokens[6]
            except IndexError:
                metacolor = None
            try:
                metanci = tokens[7]
            except IndexError:
                metanci = None
            try:
                metaumls = tokens[8]
            except IndexError:
                metaumls = None
    
            # set root node.
            prev_n = root
    
            # build nodes all the way down.
            nodes = list()
            for i in range(5):
    
                # skip empty.
                if len(tokens) < 2:
                    continue
    
                # check if empty.
                if tokens[i] == "":
                    continue
    
                # split into two.
                tmp = tokens[i].split("(")
                val = tmp[0].strip().replace('"', '').replace("'", '')
                key = tmp[1].strip().replace("(","").replace(")","").replace('"', '').replace("'", '')
    
                # build node.
                self.OT.add_node(key, {
                    'text': val,
                    'metamaintype': metamaintype,
                    'metacolor': metacolor,
                    'metanci': metanci,
                    'metaumls': metaumls
                })
                n = key
    
                # add edge.
                self.OT.add_edge(prev_n, n)
    
                # update previous node.
                prev_n = n
    
            # increment line count.
            line_cnt += 1

    
 
    def inclusion(self,cop):
        '''Function to check if cancer type exist in oncotree'''
        res = self.OT.has_node(cop)
        return res
   
 
    def homo(self, cop, coc):
        '''Function to check if two input cancer types in the same lineage'''
        while 'root' not in self.OT.predecessors(cop):
            pre_cop = self.OT.predecessors(cop)[0]
            cop = pre_cop
        while 'root' not in self.OT.predecessors(coc):
            pre_coc = self.OT.predecessors(coc)[0]
            coc = pre_coc
        return cop == coc

        
    def distance(self, cop, coc):
        '''Function to calculate distance between two input cancer types'''
        copc = 0
        cocc = 0
        if self.homo(cop, coc) == False:
            while 'root' not in self.OT.predecessors(cop):
                copc += 1
                pre_cop = self.OT.predecessors(cop)[0]
                cop = pre_cop
            while 'root' not in self.OT.predecessors(coc):
                cocc += 1
                pre_coc = self.OT.predecessors(coc)[0]
                coc = pre_coc
            total = copc + cocc + 2
            return total
        else:
            return 'cancer types are in the same lineage'
        
            

        












