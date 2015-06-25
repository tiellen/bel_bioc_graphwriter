#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python Program Description: Script for Generating a .png file containing the graph of the structure of an interaction from BioC.
If there is more than one document (document id) in the BioC file, for each document, one graph structure is generated.

Graph structures are generated with the help of graphviz (http://www.graphviz.org/) and pygraphviz (http://pygraphviz.github.io/).


"""
from optparse import OptionParser
import codecs
import sys
import os
from Bio import Entrez
import codecs
import shutil
import random
import nltk
import collections
import uuid

from bioc import BioCReader

import pygraphviz as pgv


# Prevent Encoding exceptions in Python 2.x
sys.stdout = codecs.getwriter('utf-8')(sys.__stdout__)
sys.stderr = codecs.getwriter('utf-8')(sys.__stderr__)
sys.stdin = codecs.getreader('utf-8')(sys.__stdin__)




def write_graph(bioc_graph, output_dot_dir, document_name, options=None, args=None):

    if not os.path.exists(output_dot_dir):
        os.makedirs(output_dot_dir)

    output_file_path = output_dot_dir + '/' + document_name + '.dot'

    output_dot_file = open(output_file_path, 'w')

    bioc_graph.write(output_dot_file)

    #print(bioc_graph.string())

    B = pgv.AGraph(output_file_path)

    #B = AB.acyclic(copy=True)

    picture_path = output_dot_dir + '/' + document_name + '.png'

    B.layout(prog='neato')
    B.draw(picture_path, prog='dot')
    
    os.remove(output_file_path)



def define_graph(options=None, args=None):
    one_graph = pgv.AGraph(strict=False, directed=True, ranksep='0.1')

    one_graph.node_attr['shape']='round'
    one_graph.node_attr['fixedsize']='false'
    one_graph.node_attr['fontsize']='12'
    one_graph.node_attr['style']='filled'
    one_graph.graph_attr['outputorder']='edgesfirst'
    one_graph.graph_attr['label']="BioC Graph"
    one_graph.graph_attr['ratio']='1.2'
    one_graph.edge_attr['color']='#1100FF'
    one_graph.edge_attr['style']='setlinewidth(2)'

    return one_graph


def bioc_document_to_graph(bioc_document, options=None, args=None):

        namespace_list = ['CHEBI', 'HGNC', 'MGI', 'EGID', 'GOBP', 'MESHD', 'GOCCID', '']
        pmod_arguments = ['ModificationType', 'AminoAcidCode', 'ModificationPosition', 'CodeVariant', 'Codon', 'CodeReference', 'TruncationPosition']

        bioc_graph = define_graph()

        bioc_graph.graph_attr['label']= "Graph Structure for " + bioc_document.id

        for one_passage in bioc_document.passages:

            #bioc_graph.add_nodes_from(one_passage.annotations)
            #bioc_graph.add_nodes_from(one_passage.relations)

            for one_annotation in one_passage.annotations:
                #print one_annotation
                #print one_annotation.id
                print one_annotation.infons, 'infons'
                #annotation_database_id = one_annotation.infons['database_id']
                #annotation_id_list.append(annotation_database_id)

                if 'BEL (full)' in one_annotation.infons:
                    infon_label = one_annotation.infons['BEL (full)']
                elif 'relationship' in one_annotation.infons:
                    infon_label = one_annotation.infons['relationship']
                    infon_type = one_annotation.infons['type']

                elif set(one_annotation.infons.keys()).intersection(set(namespace_list)):
                    namespace = list(set(one_annotation.infons.keys()).intersection(set(namespace_list)))[0]
                    termword = one_annotation.infons[namespace]
                    infon_label = namespace + ':' + termword
                elif set(one_annotation.infons.keys()).intersection(set(pmod_arguments)):
                    mod_arg = list(set(one_annotation.infons.keys()).intersection(set(pmod_arguments)))[0]
                    arg_code = one_annotation.infons[mod_arg]
                    infon_label = mod_arg + ':' + arg_code
                else:
                    print 'LABELNAME IS NOT DEFINED', one_annotation.infons
                    infon_label = ' '.join(one_annotation.infons.keys()) + ';' + ' '.join(one_annotation.infons.values())
                #else: annotation_label = one_annotation.infons['trigger']

                annotation_label = one_annotation.id + '\n' + infon_label

                #bioc_graph.add_node(one_annotation.id, shape='box', color='#7d9818', label=annotation_label)

                infon_type_name = one_annotation.infons['type']
                type_label = 'cluster:'+ infon_type_name + str(uuid.uuid4())
                if infon_type_name == 'relationship':
                    bioc_graph.add_node(one_annotation.id, shape='oval', color='#8B8B83', label=annotation_label)
                    bioc_graph.add_subgraph(nbunch=[one_annotation.id], name=type_label, label=infon_type_name, shape='square', color='#8B8B83', fontcolor='#8B8B83', fixedsize='false')
                else:
                    bioc_graph.add_node(one_annotation.id, shape='box', color='#7d9818', label=annotation_label)
                    bioc_graph.add_subgraph(nbunch=[one_annotation.id], name=type_label, label=infon_type_name, shape='square', color='#7d9818', fontcolor='#7d9818', fixedsize='false')

            for one_relation in one_passage.relations:
                #if the relation tag in inside the passage tag
                print one_relation
                one_relation_type = one_relation.infons['type']
                one_relation_label = one_relation.id + '\n' + one_relation_type
                bioc_graph.add_node(one_relation.id, label=one_relation_label, shape='hexagon', color='#daa623')

            for one_relation in one_passage.relations:
                for one_node in one_relation.nodes:
                    connected_node = one_node.refid
                    edge_label = one_node.role
                    bioc_graph.add_edge(one_relation.id, connected_node, label=edge_label)

        for one_relation in bioc_document.relations:
            #if the relation tag is outside the passage tag
            #first step: draw nodes
            print one_relation

            infon_label = one_relation.infons["type"]
            annotation_label = one_relation.id + '\n' + infon_label

            relation_node = bioc_graph.add_node(one_relation.id, shape='hexagon', color='#ADD8E6', label=annotation_label)
            #relation_node.node_attr.update(label=one_relation.infons["type"])


        for one_relation in bioc_document.relations:
            #second step: draw edges between nodes
            for one_node in one_relation.nodes:
                connected_node = one_node.refid
                edge_label = one_node.role
                bioc_graph.add_edge(one_relation.id, connected_node, label=edge_label)

        return bioc_graph


def process_bioc_file(bioc_file, output_dot_dir, options=None, args=None):
    '''Read BioC file with BioCReader, process each document in the BioC file seperately.'''

    bioc_reader = BioCReader(bioc_file)
    bioc_reader.read()

    for one_document in bioc_reader.collection.documents:
        doc_id = one_document.id
        print 'Document found', doc_id
        document_name = doc_id

        bioc_graph = bioc_document_to_graph(one_document)

        write_graph(bioc_graph, output_dot_dir, document_name, options=options, args=args)


def process(options=None, args=None):
    """
    Do the processing.

    The options object should be used as an argument to almost all functions.
    This gives easy access to all global parameters.
    """
    #if options.debug:
    #    print >>sys.stderr, options

    #print sys.stdin, 'test'

    print 'OPTIONS:', options

    bioc_file = args[0]
    output_dot_dir = args[1]


    graph_template = define_graph(options=options, args=args)

    process_bioc_file(bioc_file, output_dot_dir)



def main():
    """
    Invoke this module as a script
    """
    usage = "usage: %prog [options] input: BioC File (args[0]); Output Directory for the (picture) .png file."
    parser = OptionParser(version='%prog 0.99', usage=usage)

    parser.add_option('-l', '--logfile', dest='logfilename',
                      help='write log to FILE', metavar='FILE')
    parser.add_option('-q', '--quiet',
                      action='store_true', dest='quiet', default=False,
                      help='do not print status messages to stderr')
    parser.add_option('-d', '--debug',
                      action='store_true', dest='debug', default=False,
                      help='print debug information')



    (options, args) = parser.parse_args()

    if options.debug: print >> sys.stderr, '# Starting processing'

    process(options=options,args=args)




    sys.exit(0) # Everything went ok!

if __name__ == '__main__':
    main()
