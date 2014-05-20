#!/bin/env python2.7
import argparse
import os
import shutil
from bioopt.bioopt_parser import *
from routines import *

__author__ = 'Aleksej Zelezniak'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Computes nth degree CoCCoA score')
    parser.add_argument('-i', dest='infile_name', action='store',
                        help="Bioopt model file")
    parser.add_argument('-o', dest='outfile_name', action='store',
                        help="output")
    parser.add_argument('-folds', dest='exprfile_name', action='store',
                        help="Gene expression folchange pvalue file")
    parser.add_argument('-metfolds', dest='metfoldsfile_name', action='store',
                        help="Metabolite folchange pvalue file")
    parser.add_argument('-minmax', dest='mmfile_name', action='store',
                        help="Minmax file with bounds (optional)")
    parser.add_argument('-r2g', dest='r2gfile_name', action='store',
                        help="Reaction to gene(s) association file")
    parser.add_argument('-pval', dest='pval_cutoff', action='store', default=0.05, type=float,
                        help="P-value cutoff for gene expression")
    parser.add_argument('-pvalmet', dest='pvalmet_cutoff', action='store', default=1, type=float,
                        help="P-value cutoff for metabolites ")
    parser.add_argument('-degree', dest='degree', action='store', default=int(1), type=int,
                        help="Degree (int), 100 for all degrees,  default: 1")
    parser.add_argument('-met', dest='metfile_name', action='store',
                        help="List of metabolites (optional)")
    parser.add_argument('-randomize', dest='random', action='store', default=0, type=int,
                        help="Shuffles gene labels")
    parser.add_argument('-details', dest='details', action='store', default=None,
                        help="Saves genes used for calculations")
    parser.add_argument('-network', dest='network', action='store', default=None,
                        help="Shuffles gene labels")
    parser.add_argument('-forward', dest='forward', action='store_true')
    parser.add_argument('-complete', dest='complete', action='store_true')
    parser.add_argument('-direction', dest='direction', default="b", action='store')
    parser.add_argument('-verbosity', dest='verb', default=False, action='store_true')


    args = parser.parse_args()
    if not (args.infile_name and args.outfile_name):
        parser.error("Input model and output filenames are required")
    if not args.exprfile_name:
        parser.error("Fold-change pvalue expression filename is required")
    if not args.r2gfile_name:
        parser.error("Reaction to gene association filename is required")
    if not args.pval_cutoff:
        parser.error("P-value cutoff was not specified")
    if not args.degree >= 0:
        parser.error("Degree was not specified")

    infile_name = args.infile_name
    outfile_name = args.outfile_name
    exprfile_name = args.exprfile_name
    r2gfile_name = args.r2gfile_name
    pval_cutoff = args.pval_cutoff
    pvalmet_cutoff = args.pvalmet_cutoff
    degree = args.degree
    direction = args.direction
    k_times = args.random
    new_suffix = ".bak"

    metfoldsfile_name = ""
    if args.metfoldsfile_name:
        metfoldsfile_name = args.metfoldsfile_name

    mmfile_name = ""
    if args.mmfile_name:
        mmfile_name = args.mmfile_name
    metfile_name = ""
    if args.metfile_name:
        metfile_name = args.metfile_name

    detfile_name = ""
    if args.details:
        detfile_name = args.details

    netfile_name = ""
    if args.network:
        netfile_name = args.network

    warnings.simplefilter("ignore")
    model_parser = BiooptParser()
    model = model_parser.parse_file(infile_name)

    if mmfile_name:
        minmax_bounds = parse_minmax_file(mmfile_name, skip_lines=0, tolerance=1e-7)
        set_bounds(model, minmax=minmax_bounds, reverse=True)
    else:
        set_bounds(model, reverse=True)

    B = create_DiGraph(model=model, remove_rev=True, remove_nonfeas=True)

    met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0) # gets metabolites from the network
    metabolites = list(met_nodes)
    if metfile_name:
        metabolites = load_column(metfile_name)
    U = nx.projected_graph(B, met_nodes) # projecting bipartite graph to make unipartite graph of metabolites

    expression = read_expression_file(exprfile_name, sep='\s+', cutoff=pval_cutoff, folds=2, pval=1)
    reaction2gene = read_reaction2gene_file(r2gfile_name)
    results = {}
    metabolite_folds = None
    if metfoldsfile_name:
        metabolite_folds = read_metabolites_file(metfoldsfile_name, sep='\t', cutoff=pvalmet_cutoff, folds=1, pval=2)

    if args.random: #shuffling
        #random_results = []
        if os.path.exists(outfile_name):
            print "{0} exists, saving to {0}{1}".format(outfile_name, new_suffix)
            new_file = outfile_name + new_suffix
            shutil.copyfile(outfile_name, new_file)
            os.remove(outfile_name)

        for k in range(1, k_times + 1):
            if not k % 10:
                print "Past..", k
            randomized = shuffle_expression(expression)
            reaction_folds = assign_expression_nogenes(reaction2gene, randomized)
            consumption = get_consumption_scores(B, reaction_folds)
            production = get_production_scores(B, reaction_folds)
            if degree == 100:
                results = get_all_degrees_scores(U, metabolites, production, consumption)
            else:
                results = metabolites_scoring(U, metabolites, production, consumption, degree)
            print "Appending {0}...".format(k),
            save_results_append(results, outfile_name, k)
            if os.path.exists(outfile_name):
                print "OK"
    else:

        reaction_folds, reaction_genes = assign_expression(reaction2gene, expression)
        consumption = get_consumption_scores(B, reaction_folds)
        production = get_production_scores(B, reaction_folds)

        results_genes_backward = {}
        results_genes_forward = {}
        consumption_genes = {}
        production_genes = {}
        results_backward = {}
        results_forward = {}

        if detfile_name:
            consumption_genes = get_consumption_genes(B, reaction_genes)
            production_genes = get_production_genes(B, reaction_genes)

        if degree == 100 and args.forward:
            results_forward = metabolite_scoring_forward(B, metabolites, reaction_folds=reaction_folds,
                                                         metabolite_folds=metabolite_folds, degree=None)
            results_backward = metabolite_scoring_backward(B, metabolites, reaction_folds=reaction_folds,
                                                           metabolite_folds=metabolite_folds, degree=None)
            if detfile_name:
                results_genes_forward = get_all_degrees_genes_forward(U, metabolites,
                                                                      consumption_genes=consumption_genes,
                                                                      production_genes=production_genes)
                results_genes_backward = get_all_degrees_genes_backward(U, metabolites,
                                                                        consumption_genes=consumption_genes,
                                                                        production_genes=production_genes)
                results_genes = get_all_degrees_genes(U, metabolites,
                                                      consumption_genes=consumption_genes,
                                                      production_genes=production_genes)
        elif args.forward:
            results_forward = metabolite_scoring_forward(B, metabolites,
                                                         reaction_folds=reaction_folds,
                                                         metabolite_folds=metabolite_folds, degree=degree)
            results_backward = metabolite_scoring_backward(B, metabolites,
                                                           reaction_folds=reaction_folds,
                                                           metabolite_folds=metabolite_folds, degree=degree)
            if detfile_name:
                results_genes_forward = get_all_degrees_genes_forward(U, metabolites,
                                                                      consumption_genes=consumption_genes,
                                                                      production_genes=production_genes)
                results_genes_backward = get_all_degrees_genes_backward(U, metabolites,
                                                                        consumption_genes=consumption_genes,
                                                                        production_genes=production_genes)
                results_genes = get_all_degrees_genes(U, metabolites,
                                                      consumption_genes=consumption_genes,
                                                      production_genes=production_genes)
        elif degree == 100:
            results = get_all_degrees_scores(U, metabolites=metabolites,
                                             consumption=consumption, production=production)
            if detfile_name:
                results_genes_forward = get_all_degrees_genes_forward(U, metabolites,
                                                                      consumption_genes=consumption_genes,
                                                                      production_genes=production_genes)
                results_genes_backward = get_all_degrees_genes_backward(U, metabolites,
                                                                        consumption_genes=consumption_genes,
                                                                        production_genes=production_genes)
                results_genes = get_all_degrees_genes(U, metabolites,
                                                      consumption_genes=consumption_genes,
                                                      production_genes=production_genes)
        elif args.complete:
            if direction == "f":
                print "Scoring forward"
            elif direction == "b":
                print "Scoring backward"
            results = metabolite_scoring_bf(B, metabolites,
                                            reaction_folds=reaction_folds,
                                            metabolite_folds=metabolite_folds,
                                            degree=degree, dir=direction, verbose=args.verb)
        else:
            results = metabolites_scoring(U, metabolites=metabolites,
                                          consumption=consumption,
                                          production=production, degree=degree)
            if detfile_name:
                results_genes_forward = metabolites_genes_forward(U, metabolites, production_genes=production_genes,
                                                                  consumption_genes=consumption_genes, degree=degree)
                results_genes_backward = metabolites_genes(U, metabolites,
                                                           production_genes=production_genes,
                                                           consumption_genes=consumption_genes, degree=degree)
                results = get_all_degrees_genes(U, metabolites,
                                                consumption_genes=consumption_genes,
                                                production_genes=production_genes)

        if os.path.exists(outfile_name):
            print "{0} exists, saving to {0}{1}".format(outfile_name, new_suffix)
            new_file = outfile_name + new_suffix
            shutil.copyfile(outfile_name, new_file)
            os.remove(outfile_name)

        if args.forward:
            save_results(results_forward, outfile_name+".forward")
            save_results(results_backward, outfile_name+".backward")
        elif args.complete:
            if direction == "f":
                save_results(results, outfile_name+".forward")
            elif direction == "b":
                save_results(results, outfile_name+".backward")
        else:
            save_results(results, outfile_name)

        if detfile_name:
            save_results_genes(results_genes_backward, detfile_name+".backward")
            save_results_genes(results_genes_forward, detfile_name+".forward")
            save_results_genes(results_genes, detfile_name+".initial")

        if netfile_name:
            write_network_gml(B, netfile_name, reaction_folds, reaction_genes, scoring=results, degree=degree)
