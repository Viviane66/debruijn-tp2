#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Viviane"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["YViviane"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Viviane"
__email__ = "viviane.yang@etu.u-paris.fr"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.
    
    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    file = open(fastq_file, "r")
    # file_list = [i.rstrip() for i in file.read()]

    for i in file:
        yield next(file).strip()
        next(file)
        next(file)
        
        
def cut_kmer(read, k):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    """
    for line in read: 
        for pos,_ in enumerate(line[0:len(line)-k+1]):
            yield line[pos:pos+k]
    """
    read_len = len(read)
    for i in range(read_len-k+1):
        kmer=read[i:i + k]
        yield kmer
    
    
def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    dict_occur_kmers = dict()
    
    for seq in read_fastq(fastq_file):
        generator_kmers = cut_kmer(seq, kmer_size)
    
        for kmer in generator_kmers:
            if kmer in dict_occur_kmers:
                dict_occur_kmers[kmer] += 1
            else:
                dict_occur_kmers[kmer] = 1
    return dict_occur_kmers


def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    G = nx.DiGraph()
    for kmer, occurence in kmer_dict.items():
        G.add_edge(kmer[:-1], kmer[1:], weight=occurence)
    return G


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list:
        if delete_entry_node == True and delete_sink_node == True:
            # tous les noeuds
            graph.remove_nodes_from(path)
        elif delete_entry_node == True:
            # tous sauf dernier
            graph.remove_nodes_from(path[:-1])
        elif delete_sink_node == True:
            # tous sauf premier
            graph.remove_nodes_from(path[1:])
        else:
            # tous sauf premier et dernier noeud
            graph.remove_nodes_from(path[1:-1])
    return graph

        
def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    # determine best_path
    from statistics import stdev
    
    std_weight = stdev(weight_avg_list)
    if std_weight>0:
        index_path_remove = weight_avg_list.index(max(weight_avg_list))
    elif std_weight==0:
        std_length = stdev(path_length)
        
        if std_length>0:
            index_path_remove  = path_length.index(max(path_length))
        else:
            random_choice = random.randint(0, len(path_list)-1)
            index_path_remove = random_choice
            
    # remove other paths
    path_list = path_list.pop(index_path_remove)
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph
    
    
def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    # detect bubble => detects many paths
    # if many paths, calculate length and average weight
    # select the better path
    # and remove others paths 
    # return graph 
    # else, return graph 
    path_length = []
    weight_avg_list = []
        
    path_list = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    
    if len(path_list) > 0:
        for path in path_list:
            path_length.append(path)
            weight_avg_list.append(path)
        graph = select_best_path(graph, path_list, path_length, weight_avg_list)
        return graph
    else:
        return graph


def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    from itertools import combinations
    
    bubble = False
    for node in graph.nodes():
        list_predecessors = list(graph.predecessors(node))
        if len(list_predecessors) > 1:
            temp = list(combinations(list_predecessors, 2))
            for combi in list(temp):
                i, j = combi
                node_predecessor = nx.lowest_common_ancestor(graph, i, j)
                if node_predecessor != None:
                    bubble = True
                    break
        if bubble == True:
            break
    if bubble:
        graph = simplify_bubbles(solve_bubble(graph,node_predecessor, node))
    return graph
        
    
def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    start_node_delete = []
    for start_node in starting_nodes:
        list_predecessors = list(graph.predecessors(start_node)) 
        list_successors = list(graph.successors(start_node)) 
        
        if len(list_predecessors) > 0:
            start_node_delete.append(start_node) 
        elif len(list_successors) > 0:
            start_node_delete.append(start_node)
        else:
            return graph
        
    for start_node in start_node_delete:
        graph.remove_node(start_node)
    return graph
    

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    end_node_delete = []
    for end_node in ending_nodes:
        list_predecessors = list(graph.predecessors(end_node)) 
        list_successors = list(graph.successors(end_node)) 
        
        if len(list_predecessors) > 0:
            end_node_delete.append(end_node) 
        elif len(list_successors) > 0:
            end_node_delete.append(end_node)
        else:
            return graph
        
    for end_node in end_node_delete:
        graph.remove_node(end_node)
    return graph


def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    node_without_predecessors = []
    for node in graph.nodes():
        list_predecessors = list(graph.predecessors(node)) 
        if len(list_predecessors) == 0:
            node_without_predecessors.append(node)
    return node_without_predecessors
    
    
def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    node_without_successors = []
    for node in graph.nodes():
        list_successors = list(graph.successors(node)) 
        if len(list_successors) == 0:
            node_without_successors.append(node)
    return node_without_successors



def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    from itertools import product
    contig_and_length = []
    
    combi = list(product(starting_nodes, ending_nodes))
    for element in combi:
        entry_node, out_node = element
        contigs = list(nx.all_simple_paths(graph, entry_node, out_node))
        for contig in contigs:
            temp = []
            temp.append(contig)
            temp.append(len(contig))
            contig_and_length.append(temp)
    return contig_and_length
    
    
def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    with open(output_file, 'w') as file:
        for i, (contig, length) in enumerate(contigs_list, 1):
            file.write(f'>contig_{i} len={length}\n')  
            file.write(textwrap.fill(contig, width=80) + '\n')  
    return file
            


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    seqs = read_fastq(args.fastq_file)
    cut_kmer(seqs, args.kmer_size)
    build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(build_kmer_dict(args.fastq_file, args.kmer_size))
    
    # solve bubble 
    graph = simplify_bubbles(graph)
    
    # solve tips
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    
    # contigs
    #contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    #args.output_file = save_contigs(contigs_list, args.output_file)
    
    #file = open(args.output_file, "r")
    #for i in file:
    #    print(i)
    
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    if args.graphimg_file:
         draw_graph(graph, args.graphimg_file)


if __name__ == '__main__': # pragma: no cover
    main()
