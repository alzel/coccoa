import collections
from collections import Counter
import copy
import itertools
import math
import random
import re
import networkx as nx
from networkx.algorithms import bipartite
from bioopt import Bounds, Direction

__author__ = 'Aleksej Zelezniak'


def read_expression_file(path, skip_lines=0, sep="\t", id=0, folds=1, pval=2, strict=False, comment="#", cutoff=0.05,
                         filename=None, log=False):
    f = open(path, "r")
    return read_expression(f.read(), skip_lines=skip_lines, sep=sep, id=id, folds=folds, pval=pval, strict=strict,
                           comment=comment, cutoff=cutoff, filename=path, log=log)


def read_metabolites_file(path, skip_lines=0, sep="\t", id=0, folds=1, pval=2, strict=False, comment="#", cutoff=1,
                         filename=None, log=False):
    f = open(path, "r")
    return read_metabolites(f.read(), skip_lines=skip_lines, sep=sep, id=id, folds=folds, pval=pval, strict=strict,
                           comment=comment, cutoff=cutoff, filename=path, log=log)


def read_metabolites(text, skip_lines=0, sep="\t", id=0, folds=1, pval=2, strict=False, comment="#", cutoff=0.05,
                    filename=None, log=False):
    return read_expression(text=text, skip_lines=skip_lines, sep=sep, id=id, folds=folds, pval=pval, strict=strict,
                           comment=comment, cutoff=cutoff, filename=filename, log=log)


def read_expression(text, skip_lines=0, sep="\t", id=0, folds=1, pval=2, strict=False, comment="#", cutoff=0.05,
                    filename=None, log=False):
    map = {}
    error = 0
    comment_regex = re.compile("^" + comment)

    for c, cline in enumerate(text.splitlines()):
        if skip_lines > c:
            skip_lines -= 1
            continue

        cline = cline.strip()
        if not cline:
            continue
        if comment_regex.match(cline):
            continue
        parts = re.split(sep, cline)
        #parts = cline.split(sep)
        if len(parts) != 3:
            raise RuntimeError(
                "Line {0}, in {1}, must contain exactly 3 fields, but it is: {2}".format(c + 1, filename, len(parts)))
        try:
            float(parts[folds].strip())
        except:
            raise TypeError("Line {0}, in {1}, fold change must be a number".format(c + 1, filename))

        try:
            float(parts[pval].strip())
        except:
            raise TypeError("Line {0}, in {1}, p-value must be a number".format(c + 1, filename))

        fc = float(parts[folds].strip())
        p = float(parts[pval].strip())

        if p > 1 or p < 0:
            raise ValueError("Line {0}, in {1}, pvalue must be a number => 0 and <= 1".format(c + 1, filename))

        if float(p) > cutoff:
            continue

        reaction_name = parts[id].strip()

        if reaction_name in map:
            print("Line {0}, in {1}, gene is duplicated: {2}".format(c + 1, filename, reaction_name))
            error += 1

        if log and fc:
            fc = math.log(fc, 2)

        map[reaction_name] = fc

    if strict and error:
        raise RuntimeError("Found {0} critical errors".format(error))
    else:
        return map


def read_reaction2gene_file(path, skip_lines=0, sep="\t", strict=False, comment="#"):
    f = open(path, "r")
    return read_reaction2gene(f.read(), skip_lines=skip_lines, sep=sep, strict=strict, comment=comment, filename=path)


def read_reaction2gene(text, skip_lines=0, sep="\t", strict=False, comment="#", filename=None):
    map = {}
    comment_regex = re.compile("^" + comment)
    for c, cline in enumerate(text.splitlines()):
        if skip_lines > c:
            skip_lines -= 1
            continue
        cline = cline.strip()
        if not cline:
            continue

        if comment_regex.match(cline):
            continue
        parts = cline.split(sep)
        reaction_name = parts.pop(0).strip()

        if len(parts) == 0:
            print "{0} does not have gene association, in {1} ".format(reaction_name, filename)
            continue

        if reaction_name in map:
            map[reaction_name].extend([n.strip() for n in parts])
        else:
            map[reaction_name] = [n.strip() for n in parts]

    if len(map.keys()) == 0:
        raise RuntimeError("No reaction associatence were read")
    return map


def assign_expression(reaction2gene, expression, iso_pattern="\|"):
    """Assign fold change values to reaction. Averaging isoreactions. Minimum for complexes.
    """
    reaction_folds = {}
    iso_regex = re.compile(iso_pattern)
    reaction_gene = {}
    for reaction, genes in reaction2gene.iteritems():
        if len(genes) == 0:
            continue
        if re.search(iso_regex, reaction) and len(genes) > 1: #isoreactions
            tmp = []
            for gene in genes:
                if gene in expression:
                    tmp.append(float(expression[gene]))
                    if reaction in reaction_gene:
                        reaction_gene[reaction].append(gene)
                    else:
                        reaction_gene[reaction] = []
                        reaction_gene[reaction].append(gene)
                if tmp:
                    reaction_folds[reaction] = sum(tmp) / float(len(tmp)) #averaging isoreactions
        elif len(genes) > 1: #complexes
            tmp = []
            tmp_genes = [] #for storing genes in complexes in order to find which had minimal foldchange
            for gene in genes:
                if gene in expression:
                    tmp.append(float(expression[gene]))
                    tmp_genes.append(gene)
            if tmp and tmp_genes:
                reaction_folds[reaction] = min(tmp)
                gene = [gene for gene in tmp_genes if expression[gene] == reaction_folds[reaction]][0]
                if gene:
                    reaction_gene[reaction] = [gene]
                else:
                    raise RuntimeError("No gene was used for complex reaction {0}".format(reaction))
        else:
            if len(genes) != 1: #single gene case
                raise RuntimeError("Bad expression assignment for single gene reaction: ".format(reaction))
            if genes[0] in expression:
                reaction_folds[reaction] = float(expression[genes[0]])
                reaction_gene[reaction] = reaction_gene.get(genes[0], [])
                reaction_gene[reaction].append(genes[0])
    return reaction_folds, reaction_gene


def assign_expression_nogenes(reaction2gene, expression, iso_pattern="\|"):
    """Assign fold change values to reaction. Averaging isoreactions. Minimum for complexes. Does not return genes to reactions
    """
    reaction_folds = {}
    iso_regex = re.compile(iso_pattern)

    for reaction, genes in reaction2gene.iteritems():
        if len(genes) == 0:
            continue
        if re.search(iso_regex, reaction) and len(genes) > 1: #isoreactions
            tmp = []
            for gene in genes:
                if gene in expression:
                    tmp.append(float(expression[gene]))
                if tmp:
                    reaction_folds[reaction] = sum(tmp) / float(len(tmp)) #averaging isoreactions
        elif len(genes) > 1: #complexes
            tmp = []
            for gene in genes:
                if gene in expression:
                    tmp.append(float(expression[gene]))
            if tmp:
                reaction_folds[reaction] = min(tmp)
        else:
            if len(genes) != 1: #single gene case
                raise RuntimeError("Bad expression assignment for single gene reaction: ".format(reaction))
            if genes[0] in expression:
                reaction_folds[reaction] = float(expression[genes[0]])
    return reaction_folds


def save_results(results, filename):
    """Saves scores to a file.
    """
    OUT = open(filename, 'w')
    OUT.write("{0}\t{1}\t{2}\n".format("degree", "metabolite", "score"))
    for degree, mets in results.iteritems():
        for m, score in mets.iteritems():
            OUT.write("{0}\t{1}\t{2}\n".format(degree, m, score))
    OUT.close()


def save_results_genes(results, filename):
    """Saves genes to a file.
    """
    OUT = open(filename, 'w')
    OUT.write("{0}\t{1}\t{2}\n".format("degree", "metabolite", "gene"))
    for degree, mets in results.iteritems():
        for m, genes in mets.iteritems():
            for gene in genes:
                OUT.write("{0}\t{1}\t{2}\n".format(degree, m, gene))
    OUT.close()


def save_results_random(results, filename):
    """Saves scores to a file.
    """
    OUT = open(filename, 'w')
    OUT.write("{0}\t{1}\t{2}\t{3}\n".format("degree", "metabolite", "score", "trial"))
    for trial, res in enumerate(results):
        for degree, mets in res.iteritems():
            for m, score in mets.iteritems():
                OUT.write("{0}\t{1}\t{2}\t{3}\n".format(degree, m, score, trial))
    OUT.close()


def save_results_append(results, filename, trial=1):
    """Saves scores to a file.
    """
    OUT = open(filename, 'a')
    OUT.write("{0}\t{1}\t{2}\t{3}\n".format("degree", "metabolite", "score", "trial"))
    for degree, mets in results.iteritems():
        for m, score in mets.iteritems():
            OUT.write("{0}\t{1}\t{2}\t{3}\n".format(degree, m, score, trial))
    OUT.close()


def write_network_gml(B, file_name, reaction_folds, reaction_genes, scoring=None, degree=None):
    #saves network in gml format
    if degree is None:
        degree = 0
    if degree == 100 and scoring is not None:
        degree = max(scoring.keys())
    H = B.copy()
    #print reaction_folds
    met_nodes = set(n for n, d in H.nodes(data=True) if d['bipartite'] == 0)
    reaction_nodes = set(n for n, d in H.nodes(data=True) if d['bipartite'] == 1)

    for node in reaction_nodes:
        if node in reaction_genes:
            H.node[node]['label'] = "{0}:{1}".format(node, ';'.join(reaction_genes[node]))
            H.node[node]['type'] = "reaction"
        else:
            H.remove_node(node)

    for node in met_nodes:
        H.node[node]['label'] = node
        H.node[node]['type'] = "metabolite"
        if scoring and scoring[degree]:
            if node in scoring[degree]:
                H.node[node]['label'] = "{0}:{1}:{2}".format(node, degree, round(scoring[degree][node], 2))

    for edge in H.edges():
        if edge[0] in reaction_folds:
            H[edge[0]][edge[1]]['weight'] = reaction_folds[edge[0]]
        elif edge[1] in reaction_folds:
            H[edge[0]][edge[1]]['weight'] = reaction_folds[edge[1]]
    nx.write_graphml(H, file_name)


def load_column(path, key=0, skip_lines = 0, sep="\t"):
    """
    loads column, where key value denotes unique column, default is first column
    """
    array = []
    try:
        f = open(path)
        for c, cline in enumerate(f):
            if skip_lines > c:
                skip_lines -= 1
                continue

            cline = cline.strip()
            if not cline:
                continue
            parts = cline.split(sep)
            array.append(parts[key].strip())
    except IOError as e:
        print load_column.__name__,"I/O error({0}): {1}".format(e.errno, e.strerror)
        exit(-1)
    finally:
        f.close()
    return array


def parse_minmax(text, skip_lines=0, sep="\t", strict=False, comment="#", tolerance=1e-7, filename=None):
    """Parses minmax from string.
    text: (String) minmax string
    skip_lines: (Int) first n lines to skip
    sep: (String) Field delimiter
    tolerance: (Float) minimal number to consider as non-zero flux
    filename: (String) name of file, optional
    """
    map = {}
    error = 0
    comment_regex = re.compile("^" + comment)

    for c, cline in enumerate(text.splitlines()):
        if skip_lines > c:
            skip_lines -= 1
            continue

        cline = cline.strip()
        if not cline:
            continue
        if comment_regex.match(cline):
            continue

        parts = cline.split(sep)
        if len(parts) != 3:
            raise RuntimeError(
                "Line {0}, in {1}, must contain exactly 3 fields, but it is: {2}".format(c + 1, filename, len(parts)))

        reaction_name = parts.pop(0).strip()

        if reaction_name in map:
            print("Line {0}, in {1}, reaction is duplicated: {2}".format(c + 1, filename, reaction_name))
            error += 1

        lb = float(parts[0])
        ub = float(parts[1])
        if abs(lb) <= tolerance:
            lb = 0.0
        if abs(ub) <= tolerance:
            ub = 0.0
        map[reaction_name] = Bounds(lb, ub)

    if strict and error:
        raise RuntimeError("Found {0} critical errors".format(error))
    else:
        return map


def parse_minmax_file(path, skip_lines=0, sep="\t", strict=False, comment="#", tolerance=1e-7, filename=None):
    """Parses minmax from file.
    path: (String) path to minmax file
    skip_lines: (Int) first n lines to skip
    sep: (String) Field delimiter
    tolerance: (Float) minimal number to consider as non-zero flux
    filename: (String) name of file, optional
    """

    f = open(path, "r")
    return parse_minmax(f.read(), filename=path, skip_lines=skip_lines, sep=sep, strict=strict, comment=comment,
                        tolerance=tolerance)


def create_DiGraph(model, remove_rev=True, remove_nonfeas=True):
    """Creates Bipartite Directed graph from model. Nodes are created such that Reactant->Reaction->Product.
    Graph is created in the way how reactions are written in the model irrespective of Direction type.
    Reaction is considered reversible from its Bounds Direction.reversible().
    Directions must be reversed before calling this function.

    model: (bioopt.Model) model with correct bounds and reactions.
    remove_rev: (Bool) remove reversible reactions determined from bounds, default: True
    remove_nonfeas: (Bool) remove non feasible determined from bounds, default: True
    """

    G = nx.DiGraph()
    graph_reactions = []

    if remove_rev or remove_nonfeas:
        for r_i, r in enumerate(model.reactions):
            if remove_rev:
                if r.direction == Direction.reversible() and r.bounds.lb < 0 and r.bounds.ub > 0:
                    continue
            if remove_nonfeas:
                if r.bounds.lb == 0 and r.bounds.ub == 0:
                    continue
            graph_reactions.append(r)
    else:
        graph_reactions = model.reactions

    if len(graph_reactions) == 0:
        raise RuntimeError("No reactions to build graph")
    reaction_names = [r.name for r in graph_reactions]

    for r_i, r in enumerate(graph_reactions):
        r_metabolites = []
        p_metabolites = []

        for mb in r.reactants:
            if not mb.metabolite.name in reaction_names:
                r_metabolites.append(mb.metabolite.name)
            else:
                raise RuntimeError("Metabolites and reaction must have different names: {0}".format(mb.metabolite.name))
        for mb in r.products:
            if not mb.metabolite.name in reaction_names:
                p_metabolites.append(mb.metabolite.name)
            else:
                raise RuntimeError("Metabolites and reaction must have different names: {0}".format(mb.metabolite.name))
        G.add_nodes_from(r_metabolites, bipartite=0)
        G.add_nodes_from(p_metabolites, bipartite=0)
        G.add_node(r.name, bipartite=1)

        for met in r_metabolites:
            G.add_edge(met, r.name)
        for met in p_metabolites:
            G.add_edge(r.name, met)
    if bipartite.is_bipartite(G):
        return G
    else:
        raise TypeError("Cannot create bipartite graph")


def set_bounds(model, minmax=None, reverse=False):
    """Sets bounds given minmax dictionary with Bounds class values and reverses (optionally) reaction to the bounds.
    minmax: (dict) optional dictionaty with bounds, otherwise will use as provided in the model
    reverse: (bool) reverses reaction in the model if set to True and if lower bounds < 0 and upper <= 0
    """
    if minmax is None and reverse:
        for r in model.reactions:
            if reverse and r.direction == Direction.reversible() and r.bounds.lb < 0 and r.bounds.ub <= 0:
                r.reverse()
        return model
    elif minmax is None:
        print "Nothing done to bounds"
        return model
    else:
        if not isinstance(minmax, dict):
            raise TypeError("minmax must be dictionary")

        if len(minmax.keys()) == 0:
            raise RuntimeError("minmax should contain at least 1 reaction")

        for val in minmax.values():
            if not isinstance(val, Bounds):
                raise TypeError("minmax values must be bounds")

        for r in model.find_reactions(minmax.keys()):
            lb = minmax[r.name].lb
            ub = minmax[r.name].ub
            r.bounds.lb = lb
            r.bounds.ub = ub
            if reverse and r.direction == Direction.reversible() and r.bounds.lb < 0 and r.bounds.ub <= 0:
                r.reverse()
        return model


def get_output_metabolites(B, n):
    return [l for m in B.successors(n) for l in B.successors(m)]


def get_input_metabolites(B, n):
    return [l for m in B.predecessors(n) for l in B.predecessors(m)]


def get_consumption_scores(B, reaction_folds):
    """For each metabolite in B returns dictionary average reaction fold-changes from consumption side
    """
    scores = {}
    met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
    for metabolite in met_nodes:
        tmp_score = []
        for reaction in B.successors(metabolite):
            if reaction in reaction_folds:
                tmp_score.append(reaction_folds[reaction])
        if tmp_score:
            scores[metabolite] = float(sum(tmp_score)) / float(len(tmp_score))
        else:
            scores[metabolite] = float(0)
    if not scores:
        raise RuntimeError("None of consumption scores were computed")
    return scores


def get_production_scores(B, reaction_folds):
    """For each metabolite in B returns dictionary average reaction fold-changes from production side
    """
    scores = {}
    met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
    for metabolite in met_nodes:
        tmp_score = []
        for reaction in B.predecessors(metabolite):
            if reaction in reaction_folds:
                tmp_score.append(reaction_folds[reaction])
        if tmp_score:
            scores[metabolite] = float(sum(tmp_score)) / float(len(tmp_score))
        else:
            scores[metabolite] = float(0)
    if not scores:
        raise RuntimeError("None of production scores were computed")
    return scores


def get_consumption_genes(B, reaction_genes):
    """For each metabolite in B returns dictionary of arrays with consumption genes of each metabolite
    """
    genes = {}
    met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
    for metabolite in met_nodes:
        tmp_genes = []
        for reaction in B.successors(metabolite):
            if reaction in reaction_genes:
                tmp_genes.extend(reaction_genes[reaction])
        if tmp_genes:
            genes[metabolite] = list(set(tmp_genes))
    return genes


def get_production_genes(B, reaction_genes):
    """For each metabolite in B returns dictionary of arrays with production genes of each metabolite
    """
    genes = {}
    met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
    for metabolite in met_nodes:
        tmp_genes = []
        for reaction in B.predecessors(metabolite):
            if reaction in reaction_genes:
                tmp_genes.extend(reaction_genes[reaction])
        if tmp_genes:
            genes[metabolite] = list(set(tmp_genes))
    return genes


def get_production_scores2(B, consumption):
    """For each metabolite in B returns dictionary average reaction fold-changes from production side
    """
    scores = {}
    met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
    U = nx.projected_graph(B, met_nodes)
    Urev = U.reverse()

    for metabolite in met_nodes:
        path = nx.single_source_dijkstra_path_length(Urev, metabolite, cutoff=1)
        n_score = []
        for neighbour in path:
            if path[neighbour] == 1:
                if neighbour in consumption and consumption[neighbour] != 0:
                    n_score.append(consumption[neighbour])
        if n_score:
            scores[metabolite] = float(sum(n_score)) / float(len(n_score))
        else:
            scores[metabolite] = float(0)

    if not scores:
        raise RuntimeError("None of production scores were computed")
    return scores


def metabolites_scoring(U, metabolites, consumption, production, degree=None):
    """For a given list of metabolites and degree calculates CoCCoA scores. By default calculates maximum possible degree.
    U: (DiGraph) graph of metabolites.
    metabolites: (list) list of metabolites.
    consumption/production: (dict) average consumption/production scores for all metabolites
    degree: (int) CoCCoA degree
    """
    if not isinstance(metabolites, list):
        raise TypeError("metabolites must be list")
    if not U.is_directed():
        raise TypeError("Graph must be directed")

    effective_degree = 0
    if degree >= 1:
        effective_degree = degree - 1 #effective degree since it is actually distance cuttoff
    if degree is None:
        effective_degree = None

    found_metabolites = list(set(U.nodes()) & set(metabolites))
    not_found = list(set(metabolites) - set(U.nodes()))

    if not_found:
        for m in not_found:
            print "No such metabolite {0}".format(m)

    if len(found_metabolites) == 0:
        raise TypeError("No metabolites found in the network")

    scores = collections.defaultdict(dict) #final results
    Urev = U.reverse() #reversing graph to get correct distances from source

    for metabolite in found_metabolites:
        path = nx.single_source_dijkstra_path_length(Urev, metabolite,
                                                     cutoff=effective_degree) #returns nodes and shortest path lengths until nodes, cutoff is maximal distance for search
        distances = {}

        for k, v in path.iteritems(): #reversing path lengths to get nodes which are at the distance from source metabolite distance[degree][metabolite]
            if v == 0 and degree == 0: #for zero degree
                if 0 in distances:
                    distances[0].append(k)
                else:
                    distances[0] = [k]
                break
            if v == 0: #for all the rest degrees
                if 0 in distances:
                    distances[0].append(k)
                else:
                    distances[0] = [k]
                if 1 in distances:
                    distances[1].append(k)
                else:
                    distances[1] = [k]
            elif v >= 1:
                n = v + 1 #
                if n in distances:
                    distances[n].append(k)
                else:
                    distances[n] = [k]

        score = float(0)
        for distance, mets in distances.iteritems():
            dist_score = []
            if distance == 0: #special case
                cons_m = float(0)
                if metabolite in consumption:
                    cons_m = consumption[metabolite]
                if cons_m != 0:
                    dist_score.append(-cons_m)
            elif distance == 1: #special case
                #cons_m = float(0)
                prod_m = float(0)
                # if metabolite in consumption:
                #     cons_m = consumption[metabolite]
                if metabolite in production:
                    prod_m = production[metabolite]
                if prod_m != 0:
                    dist_score.append(prod_m)
            else:
                for m in mets:
                    cons_m = float(0)
                    prod_m = float(0)
                    if m in consumption:
                        cons_m = consumption[m]
                    if m in production:
                        prod_m = production[m]
                    if len(mets) > 1 and cons_m != 0 or len(mets) > 1 and prod_m != 0: # average (-cons + prod) per degree in order to avoid 0 score for averaging
                        dist_score.append(-cons_m + prod_m)
                    else:
                        dist_score.append(-cons_m + prod_m)
            if dist_score:
                score += sum(dist_score) / float(len(dist_score))
        scores[max(distances.keys())][metabolite] = score
    return scores


def metabolite_scoring_forward(B, metabolites, reaction_folds, metabolite_folds = None,  degree=None):
    """For a given list of metabolites and degree returns genes for CoCCoA scores. By default calculates maximum possible degree.
    U: (DiGraph) graph of metabolites.
    metabolites: (list) list of metabolites.
    degree: (int) CoCCoA degree
    """
    if not isinstance(metabolites, list):
        raise TypeError("metabolites must be list")
    if not bipartite.is_bipartite(B):
        raise TypeError("Graph must be bipartite")

    met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
    U = nx.projected_graph(B, met_nodes) # projecting bipartite graph to make unipartite graph of metabolites

    found_metabolites = list(set(U.nodes()) & set(metabolites))
    not_found = list(set(metabolites) - set(U.nodes()))

    if not_found:
        for m in not_found:
            print "No such metabolite {0}".format(m)

    if len(found_metabolites) == 0:
        raise TypeError("No metabolites found in the network")

    scores = collections.defaultdict(dict) #final results
    average_metabolite = 0

    if metabolite_folds and len(metabolite_folds.values()):
        average_metabolite = sum(metabolite_folds.values())/float(len(metabolite_folds.values()))


    for metabolite in found_metabolites:
        path_lengths = nx.single_source_dijkstra_path_length(U, metabolite,
                                                             cutoff=degree) #returns nodes and shortest path lengths until nodes, cutoff is maximal distance for search
        #print path_lengths
        distances = {}
        for k, v in path_lengths.iteritems(): #reversing path lengths to get nodes which are at the distance from source metabolite distance[distance][metabolite]
            if v in distances:
                distances[v].append(k)
            else:
                distances[v] = [k]
        distances = collections.OrderedDict(sorted(distances.items()))
        path_members = {}
        metabolite_paths = set()
        #print distances
        reaction_edges = set() # contains reactions used in paths
        path_hash = collections.defaultdict(dict) #stores hash[distance][source][target] =  reaction

        pathways = {} #stores shortest paths at distance

        for distance, mets in distances.iteritems():
            #print "Distance metabolites:", distance, mets
            pathways[distance] = []
            for m in mets:
                shortest_paths = [p for p in nx.all_shortest_paths(U, source=metabolite, target=m)]
                for path in shortest_paths:
                    pathways[distance].append(path)
                    source = 0
                    target = source + 1
                    #print path
                    while target < len(path): #collecting all reactions connecting metabolites for each path
                        #print path[source], path[target]
                        common_edges = set()
                        s_node = path[source]
                        t_node = path[target]
                        source_cons = B.successors(s_node)
                        target_prod = B.predecessors(t_node)
                        metabolite_paths.add(s_node)
                        metabolite_paths.add(t_node)

                        #print "source node, reaction", s_node, source_cons
                        #print "target node, reaction", t_node, target_prod
                        common_edges = set(target_prod) & set(source_cons) # common edge
                        #print "Common", common_edges
                        #print "######################"

                        if distance in path_members:
                            path_members[distance].update(common_edges)
                        else:
                            path_members[distance] = common_edges

                        if len(common_edges) < 1:
                            raise RuntimeError("Something REALLY wrong!")
                        #print "Source {0}, target {1}".format(s_node, t_node)
                        #path_hash[s_node][t_node] = common_edges.copy()
                        path_hash[s_node][t_node] = common_edges.copy()

                        # for common_reaction in common_edges:
                        #     path_hash[s_node][t_node].add(common_reaction)
                        reaction_edges.update(common_edges)
                        source = target
                        target += 1
        #print pathways

        for distance_global in distances.keys():
            pathways_distance = pathways.copy()
            max_distance = distance_global
            current = 0
            next = current + 1
            while current < max_distance: #filtering pathways , removing pathway if those are subpathways
                for i, current_path in enumerate(pathways_distance[current]):
                    c_path = set(current_path)
                    for next_path in pathways_distance[next]:
                        n_path = set(next_path)
                        if c_path.issubset(n_path):
                            pathways_distance[current].pop(i)
                            break
                current = next
                next += 1
            for i in pathways_distance.keys():
                if i > distance_global:
                    del pathways_distance[i]

            #print pathways
            #print path_hash
            dist_score = []
            accounted = {}
            pathways_distance = collections.OrderedDict(sorted(pathways_distance.items()))
            #print "distance, pathway", distance_global, pathways_distance
            for distance, paths in pathways_distance.iteritems(): #notice paths are in reverse order
                if distance == 0:
                    continue
                if len(paths) == 0:
                    continue
                for path in paths:
                    denominator = 1
                    len_path = float(len(path)) - 1.0

                    if not len_path:
                        raise RuntimeError("Something REALLY wrong, pathway should have at least one node!")
                    if not denominator:
                        raise RuntimeError("Something REALLY wrong, node metabolite seems to be not connected!")

                    #print "starting scoring in path",  distance, path
                    source = 0
                    target = source + 1
                    while target < len(path):
                        s_node = path[source]
                        t_node = path[target]
                        #print "Source {0}, target {1}".format(s_node, t_node)

                        if s_node in path_hash and path_hash[s_node][t_node]:
                            node_production_degree = float(U.out_degree(s_node))
                            #node_production_degree = float(B.out_degree(s_node))

                            #tmp_degree = 0
                            #for r in B.successors(s_node):
                            #    if r in reaction_folds:
                            #        tmp_degree += 1
                            #if tmp_degree:
                            #    node_production_degree = float(tmp_degree)

                            if node_production_degree:
                                denominator *= node_production_degree

                            for prod_reaction in path_hash[s_node][t_node]:
                                #print "Production reaction of s_node", s_node, prod_reaction
                                if prod_reaction in accounted: #checking if reaction was already counted
                                   continue
                                else:
                                    # prod_reaction
                                    accounted[prod_reaction] = 1

                                reaction_score = 0 #score of each flux
                                if prod_reaction in reaction_folds:
                                    reaction_score = 1.0/(float(len_path)*denominator)*reaction_folds[prod_reaction]/len(path_hash[s_node][t_node])
                                    #print prod_reaction, reaction_score, denominator, (float(len_path)*denominator)

                                dist_score.append(reaction_score)

                                metabolite_score = 0 #adding metabolite fold changes information, if available
                                if metabolite_folds and s_node != metabolite and s_node in metabolite_folds:
                                    metabolite_score = 1.0/(float(len_path)*denominator)*metabolite_folds[s_node]/len(path_hash[s_node][t_node])
                                elif s_node != metabolite and average_metabolite:
                                    metabolite_score = 1.0/(float(len_path)*denominator)*average_metabolite/len(path_hash[s_node][t_node])
                                dist_score.append(metabolite_score)

                                if target + 1 != len(path): #checking if we are not at the edge of graph
                                    minus_reactions = set() #minusing consuming/producing reactions depending whether backward or forward
                                    node_consumers = B.predecessors(t_node) # consumption of target node, notice graph is reversed
                                    minus_reactions = set(node_consumers) - reaction_edges # minus reaction which are consumption of source not but are not used in shortest paths
                                    metabolite_minus = set(U.predecessors(t_node)) - metabolite_paths #getting metabolites which are precursors to t_node in this path, but are absent from pathway


                                    #print minus_reactions
                                    if minus_reactions and max_distance > 1: #minus consumption
                                        node_consumption_degree = float(len(minus_reactions))
                                        denominator_minus = denominator*node_consumption_degree
                                        for reaction in minus_reactions:
                                            if reaction in reaction_folds:
                                                minus_reaction_score = 1.0/(float(len_path)*denominator_minus)*reaction_folds[reaction]
                                                #print reaction, minus_reaction_score
                                                dist_score.append(-minus_reaction_score)

                                            for met in metabolite_minus:
                                                metabolite_score = 0
                                                if metabolite_folds and met in metabolite_folds:
                                                    metabolite_score = 1.0/(float(len_path)*denominator_minus)*metabolite_folds[met]/len(path_hash[s_node][t_node])
                                                elif average_metabolite:
                                                    metabolite_score = 1.0/(float(len_path)*denominator_minus)*average_metabolite/len(path_hash[s_node][t_node])
                                                dist_score.append(-metabolite_score)


                            source = target
                            target += 1

            scores[distance_global][metabolite] = sum(dist_score)
    return scores


def find_sublist(sub, bigger):
    if not bigger:
        return -1
    if not sub:
        return 0
    first, rest = sub[0], sub[1:]
    pos = 0
    try:
        while True:
            pos = bigger.index(first, pos) + 1
            if not rest or bigger[pos:pos+len(rest)] == rest:
                return pos -1
    except ValueError:
        return -1


def is_sublist(sub, bigger):
    if not sub:
        return False
    else:
        return find_sublist(sub, bigger) >= 0


def _filter_pathways(pathways, max_distance):
    """For a given unique pathways per distance (pathways), removes pathways which are subsets of pathways at distance and at the next
    distance until max_distance.
    Returns filtered dictionary until max_distance.
    :param pathways:
    :param max_distance:
    pathways_reactions: (dict[distance][[path1][path2][path3]) dictionary of pathways
    max_distance: (int) maximum distance until to filter
    """
    pathways_distance = copy.deepcopy(pathways) #pathways at given distance
    if max_distance > max(pathways_distance.keys()):
        max_distance = max(pathways_distance.keys())
    #removing subsets of pathways at current degree.
    for d, paths in pathways_distance.iteritems():
        if not paths:
            continue
        for i in xrange(len(paths)-1,-1,-1):
            for j in xrange(i-1, -1, -1):
                if paths[j][0] != paths[i][0]:
                    continue
                bigger = []
                sub = []
                if len(paths[j]) == len(paths[i]):
                    continue
                elif len(paths[j]) < len(paths[i]):
                    bigger = paths[i]
                    sub = paths[j]
                else:
                    bigger = paths[j]
                    sub = paths[i]
                #print paths[i], paths[j]
                if is_sublist(sub, bigger):
                    #print "removing", sub
                    pathways_distance[d].remove(sub)
                    break
    #pathways_new = copy.deepcopy(pathways_distance)
    current = 0
    next = current + 1
    while current < max_distance: #filtering pathways, removing pathway if those are subsets of +1 more distant pathways
        for i in xrange(len(pathways_distance[current])-1,-1,-1):
            current_path = pathways_distance[current][i]
            sub = current_path
            for next_path in pathways_distance[next]:
                if current_path[0] != next_path[0]:
                    continue
                if current_path == next_path:
                    pathways_distance[current].remove(sub)
                    break
                bigger = next_path
                if is_sublist(sub, bigger):
                    pathways_distance[current].remove(sub)
                    break
        current = next
        next += 1

    for i in pathways_distance.keys(): #removing rest pathways those that are larger than max_distance
        if i > max_distance:
            del pathways_distance[i]
    return pathways_distance


def _filter_pathways_light(pathways, max_distance):
    pathways_distance = copy.deepcopy(pathways)
    current = 0
    next = current + 1
    while current < max_distance: #filtering pathways, removing pathway if those are subsets of +1 more distant pathways
        for i in xrange(len(pathways_distance[current])-1,-1,-1):
            current_path = pathways_distance[current][i]
            sub = current_path
            for next_path in pathways_distance[next]:
                if current_path[0] != next_path[0]:
                    continue
                if current_path == next_path:
                    pathways_distance[current].remove(sub)
                    break
                bigger = next_path
                if is_sublist(sub, bigger):
                    pathways_distance[current].remove(sub)
                    break
        current = next
        next += 1

    for i in pathways_distance.keys(): #removing rest pathways those that are larger than max_distance
        if i > max_distance:
            del pathways_distance[i]
    return pathways_distance


def metabolite_scoring_bf(B, metabolites, reaction_folds, metabolite_folds=None, degree=None, dir="b", verbose = True):
    """For a given list of metabolites and degree returns genes for CoCCoA scores. By default calculates maximum possible degree.
    :param B:
    :param metabolites:
    :param reaction_folds:
    :param metabolite_folds:
    :param degree:
    U: (DiGraph) graph of metabolites.
    metabolites: (list) list of metabolites.
    degree: (int) CoCCoA degree
    cons1: (bool) if first degree consumtion shuold be substracted
    """

    if not isinstance(metabolites, list):
        raise ValueError("metabolites must be list")
    if not bipartite.is_bipartite(B):
        raise TypeError("Graph must be bipartite")
    if dir not in ["f", "b"]:
        raise ValueError("Argument dir must be 'b' or 'f'")
    met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
    react_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 1)
    U = nx.projected_graph(B, met_nodes) # projecting bipartite graph to make unipartite graph of metabolites
    found_metabolites = list(set(U.nodes()) & set(metabolites))
    not_found = list(set(metabolites) - set(U.nodes()))

    if not_found:
        for m in not_found:
            print "No such metabolite {0}".format(m)

    if len(found_metabolites) == 0:
        raise ValueError("No metabolites found in the network")

    scores = collections.defaultdict(dict) #final results

    G = None
    if dir == "f":
        G = U.copy()
    else:
        G = U.reverse() # reversing graph to get correct distances from source


    average_metabolite = 0
    if metabolite_folds and len(metabolite_folds.values()):
        average_metabolite = sum(metabolite_folds.values())/float(len(metabolite_folds.values()))

    for metabolite in found_metabolites:
        path_lengths = nx.single_source_dijkstra_path_length(G, metabolite, cutoff=degree)  # returns nodes and shortest path lengths until nodes, cutoff is maximal distance for search
        distances = {}
        for k, v in path_lengths.iteritems():  # reversing path lengths to get nodes which are at the distance from source metabolite distance[distance][metabolite]
            if v in distances:
                distances[v].append(k)
            else:
                distances[v] = [k]

        distances = collections.OrderedDict(sorted(distances.items()))
        path_members = {}
        metabolite_paths = set()
        reaction_edges = set()  # contains reactions used in paths
        metabolite_nodes = set()  # contains metabolites used in pathways
        path_hash = collections.defaultdict(dict)  # stores hash[distance][source][target] =  reaction
        #effective_paths = collections.defaultdict(dict)

        pathways = {}  # stores shortest paths at distance
        pathways_reactions = {}  # stores shortest paths at distance in reactions in form of dictionary[distance] = [[set(r1),set(r2,r3)], [set(r1),set(r4,r5)]]
        total = 0
        for distance, mets in distances.iteritems():
            #print "Distance metabolites:", distance, mets
            pathways[distance] = []
            pathways_reactions[distance] = []

            for m in mets:
                #if verbose:
                #    print "all simple paths at from {0}, to {1}, at distance {2}...".format(metabolite, m, distance+1)
                shortest_paths = [p for p in nx.all_simple_paths(G, source=metabolite, target=m, cutoff=distance+1)]
                total += len(shortest_paths)
                for path in shortest_paths:
                    pathways[distance].append(path)
                    source = 0
                    target = source + 1
                    reactions_path = []  # stores assembled path of reactions
                    r_paths = []
                    while target < len(path):  # collecting all reactions connecting metabolites for each path
                        s_node = path[source]  # current node
                        t_node = path[target]  # next node
                        metabolite_nodes.update(s_node)
                        metabolite_nodes.update(t_node)
                        source_cons = []
                        target_prod = []
                        if dir == "f":
                            source_cons = B.successors(s_node)
                            target_prod = B.predecessors(t_node)
                        else:
                            source_cons = B.predecessors(s_node)
                            target_prod = B.successors(t_node)

                        common_edges = set(target_prod) & set(source_cons)  # common edge

                        reaction_edges.update(common_edges)
                        if len(common_edges) < 1:
                            raise RuntimeError("Something REALLY wrong!")
                        path_hash[s_node][t_node] = common_edges.copy()
                        if len(common_edges) == 1:
                            if not r_paths:
                                r_paths.append([])
                            for i in range(0, len(r_paths)):
                                r_paths[i].append(list(common_edges)[0])
                        else:  # makes number of copies of current path and extends it to include multiple common_edges
                            if not r_paths:
                                for i in list(common_edges):
                                    r_paths.append([i])
                            else:
                                tmp = []
                                for p in r_paths:
                                    for i in itertools.product(p, list(common_edges)): #
                                        tmp.append(list(i))
                                r_paths.extend(tmp)
                        source = target
                        target += 1
                    pathways_reactions[distance].extend(r_paths)
        if verbose:
            print "{0} total pathways: {1} at distance {2}".format(metabolite, total, degree)


        uniq_pathways_reactions = collections.defaultdict(list)  # for storing unique reaction pathways per distance
        for d, list_p in pathways_reactions.iteritems():  # removing duplicated pathways e.g. resulted from m1->m2 + m3
            uniq_p_set = set(map(tuple, list_p))
            uniq_pathways_reactions[d] = map(list, uniq_p_set)

        #print "###Path_hash", path_hash
        #print distances.keys()
        #print "Pathways", pathways
        #print "before filtering", uniq_pathways_reactions
        for distance_global in distances.keys():
            if distance_global == 0:
                continue
            #print "For degree", distance_global
            pathways_distance = _filter_pathways(pathways=uniq_pathways_reactions, max_distance=distance_global)
            #print "After filtering", pathways_distance
            total_pathways = 0
            reaction_l_path = {}  # stores shortest path where reaction was present
            for d, list_p in pathways_distance.iteritems():
                #print list_p
                total_pathways += len(list_p)
                for path in list_p:
                    for r in path:
                        if r in reaction_l_path:
                            if reaction_l_path[r] > len(path):
                                reaction_l_path[r] = len(path)
                        else:
                            reaction_l_path[r] = len(path)
            if not total_pathways:
                raise RuntimeError("Something REALLY wrong: no pathways found at distance {0}!".format(distance_global))

            reaction_counts = Counter()  # counts of every reaction per distance
            reaction_weights = Counter()  # reaction final weights

            for d, list_p in pathways_distance.iteritems():
                if d == 0 or len(list_p) == 0:
                    continue
                for path in list_p:
                    for r in path:
                        if r in reaction_l_path:
                            reaction_counts[r] += 1
                            reaction_weights[r] += 1/float(reaction_l_path[r])*1/float(total_pathways)

            #print reaction_counts
            #print reaction_weights
            #print reaction_l_path

            dist_score = []
            accounted = {}
            accounted_minus = {}
            accounted_metabolites = {}
            accounted_metabolites_minus = {}
            first_forward_neighbours = B.successors(metabolite)
            for distance, paths in pathways_distance.iteritems():
                if distance == 0:
                    continue
                if len(paths) == 0:
                    continue
                for path in paths:
                    for i in range(0, len(path)):
                        reaction = path[i]
                        reaction_score = 0  # scoring each reaction
                        if reaction in reaction_folds and reaction not in accounted:
                            reaction_score = reaction_weights[reaction]*reaction_folds[reaction]
                            #print "Reaction: {0}, Weight: {1}, Folds: {2}, Score: {3}".format(reaction, reaction_weights[reaction], reaction_folds[reaction], reaction_score)
                            accounted[reaction] = 1
                        dist_score.append(reaction_score)

                        #adding metabolite fold changes information, if available
                        metabolite_score = 0
                        if metabolite_folds:
                            met_comp = B.predecessors(reaction)
                            for met in met_comp: # if there are more than one metabolite per reaction take average
                                if dir == "f" and reaction in first_forward_neighbours: #skipping first neighbours in forward direction
                                    #print "not accounting", reaction, met
                                    break
                                #print reaction, met
                                if met in metabolite_folds and (reaction, met) not in accounted_metabolites:
                                    metabolite_score += reaction_weights[reaction]*metabolite_folds[met]/float(len(met_comp))
                                    accounted_metabolites[reaction, met] = 1
                                elif average_metabolite and (reaction, met) not in accounted_metabolites:
                                    metabolite_score += reaction_weights[reaction]*average_metabolite/float(len(met_comp))
                                    accounted_metabolites[reaction, met] = 1
                                #print metabolite_score
                        dist_score.append(metabolite_score)

                        #minusing consumer reactions of target metabolite
                        if i < len(path) - 1:  # makes sure that it is not the last flux in a pathway
                            #print i, reaction
                            mets_infront = []
                            if dir == "f":
                                mets_infront = B.successors(reaction) # metabolites in front of current flux, notice scoring graph is reversed
                            else:
                                mets_infront = B.predecessors(reaction)
                            #mets_infront = B.predecessors(reaction)
                            #print reaction, mets_infront
                            for met in mets_infront:
                                tmp_reactions = []
                                if dir == "f":
                                    tmp_reactions = B.predecessors(met) # metabolites in front of current flux, notice scoring graph is reversed
                                else:
                                    tmp_reactions = B.successors(met)
                                #tmp_reactions = B.successors(met)
                                minus_reactions = set(tmp_reactions) - reaction_edges  # making sure reaction is not a part of paths
                                minus_metabolite_score = 0
                                minus_reaction_score = 0
                             #   print minus_reactions
                                for minus_r in minus_reactions:
                                    if minus_r in reaction_folds and (reaction, minus_r) not in accounted_minus:
                                        minus_reaction_score -= reaction_folds[minus_r]*reaction_weights[reaction]
                                        accounted_minus[reaction, minus_r] = 1
                                    if metabolite_folds:
                                        tmp_metabolites = B.predecessors(minus_r)
                                        minus_metabolites = set(tmp_metabolites)
                                        for minus_m in minus_metabolites:
                                            if dir == "f" and minus_r in first_forward_neighbours: #skipping first neighbours in forward direction
                                                break
                                            if minus_m in metabolite_folds and (reaction, minus_r, minus_m) not in accounted_metabolites_minus:
                                                minus_metabolite_score -= reaction_weights[reaction]*metabolite_folds[minus_m]/float(len(minus_metabolites))
                                                accounted_metabolites_minus[reaction, minus_r, minus_m] = 1
                                            elif average_metabolite and (reaction, minus_r, minus_m) not in accounted_metabolites_minus:
                                                minus_metabolite_score -= reaction_weights[reaction]*average_metabolite/float(len(minus_metabolites))
                                                accounted_metabolites_minus[reaction, minus_r, minus_m] = 1
                                dist_score.append(minus_reaction_score)
                                dist_score.append(minus_metabolite_score)
            scores[distance_global][metabolite] = sum(dist_score)
    return scores


def metabolite_scoring_backward(B, metabolites, reaction_folds, metabolite_folds=None, degree=None):
    """For a given list of metabolites and degree returns genes for CoCCoA scores. By default calculates maximum possible degree.
    U: (DiGraph) graph of metabolites.
    metabolites: (list) list of metabolites.
    degree: (int) CoCCoA degree
    """

    if not isinstance(metabolites, list):
        raise TypeError("metabolites must be list")
    if not bipartite.is_bipartite(B):
        raise TypeError("Graph must be bipartite")

    met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
    U = nx.projected_graph(B, met_nodes) # projecting bipartite graph to make unipartite graph of metabolites

    found_metabolites = list(set(U.nodes()) & set(metabolites))
    not_found = list(set(metabolites) - set(U.nodes()))

    if not_found:
        for m in not_found:
            print "No such metabolite {0}".format(m)

    if len(found_metabolites) == 0:
        raise TypeError("No metabolites found in the network")

    scores = collections.defaultdict(dict) #final results
    Urev = U.reverse() #reversing graph to get correct distances from source

    average_metabolite = 0

    if metabolite_folds and len(metabolite_folds.values()):
        average_metabolite = sum(metabolite_folds.values())/float(len(metabolite_folds.values()))

    for metabolite in found_metabolites:
        path_lengths = nx.single_source_dijkstra_path_length(Urev, metabolite,
                                                             cutoff=degree) #returns nodes and shortest path lengths until nodes, cutoff is maximal distance for search
        distances = {}
        for k, v in path_lengths.iteritems(): #reversing path lengths to get nodes which are at the distance from source metabolite distance[distance][metabolite]
            if v in distances:
                distances[v].append(k)
            else:
                distances[v] = [k]
        distances = collections.OrderedDict(sorted(distances.items()))
        path_members = {}

        #print "Source metabolite:", metabolite, "Distances:", distances
        metabolite_paths = set()
        reaction_edges = set() # contains reactions used in paths
        path_hash = collections.defaultdict(dict) #stores hash[distance][source][target] =  reaction
        #effective_paths = collections.defaultdict(dict)
        #print distances
        pathways = {} #stores shortest paths at distance

        for distance, mets in distances.iteritems():
            #print "Distance metabolites:", distance, mets
            pathways[distance] = []
            for m in mets:
                #tmp_path_hash = collections.defaultdict(dict)
                shortest_paths = [p for p in nx.all_shortest_paths(Urev, source=metabolite, target=m)]

                for path in shortest_paths:
                    pathways[distance].append(path)
                    source = 0
                    target = source + 1
                    while target < len(path): #collecting all reactions connecting metabolites for each path
                        #print path[source], path[target]
                        common_edges = set()
                        s_node = path[source]
                        t_node = path[target]


                        source_cons = B.predecessors(s_node)
                        target_prod = B.successors(t_node)

                        #print "source node, reaction", s_node, source_cons
                        #print "target node, reaction", t_node, target_prod
                        common_edges = set(target_prod) & set(source_cons) # common edge
                        if len(common_edges) < 1:
                            raise RuntimeError("Something REALLY wrong!")
                        #print "Source {0}, target {1}".format(s_node, t_node)
                        path_hash[s_node][t_node] = common_edges.copy()

                        # for common_reaction in common_edges:
                        #     path_hash[s_node][t_node].add(common_reaction)
                        reaction_edges.update(common_edges)
                        source = target
                        target += 1
        #print "###Path_hash", path_hash
        print distances.keys()
        for distance_global in distances.keys():
            print "For degree", distance_global

            pathways_distance = _filter_pathways_light(pathways=pathways, max_distance=distance_global)
            dist_score = []
            accounted = {}
            accounted_metabolites = collections.defaultdict(dict)
            pathways_distance = collections.OrderedDict(sorted(pathways_distance.items()))

            #print "Global distance {0}, pathway{1}".format(distance_global, pathways_distance)
            for distance, paths in pathways_distance.iteritems(): #notice paths are in reverse order
                if distance == 0:
                    continue
                if len(paths) == 0:
                    continue
                for path in paths:
                    denominator = 1
                    len_path = float(len(path)) - 1.0 #because we are interested in number of edges
                    if not len_path:
                        raise RuntimeError("Something REALLY wrong, pathway should have at least one node!")
                    if not denominator:
                        raise RuntimeError("Something REALLY wrong, node metabolite seems to be not connected!")

                    print "starting scoring in path",  distance, path
                    source = 0
                    target = source + 1
                    while target < len(path):
                        s_node = path[source]
                        t_node = path[target]
                        #print "Source {0}, target {1}".format(s_node, t_node)

                        if s_node in path_hash and path_hash[s_node][t_node]:

                            node_split_degree = float(B.in_degree(s_node))
                            out_reactions = B.successors(s_node)
                            used_edges = float(len(set(out_reactions) & common_edges))  #
                            if used_edges > 1:
                                node_split_degree /= used_edges

                            if node_split_degree:
                                denominator *= node_split_degree

                            for prod_reaction in path_hash[s_node][t_node]:
                               # print "Production reaction of", s_node, " is ", prod_reaction
                                if prod_reaction in accounted: #checking if reaction was already counted and its metabolites
                                    if prod_reaction in accounted_metabolites and metabolite_folds and t_node not in accounted_metabolites[prod_reaction]:
                                        metabolite_score = 0
                                        if metabolite_folds and t_node in metabolite_folds:
                                            metabolite_score = 1.0/(float(len_path)*denominator)*metabolite_folds[t_node]/len(path_hash[s_node][t_node])
                                        elif average_metabolite:
                                            metabolite_score = 1.0/(float(len_path)*denominator)*average_metabolite/len(path_hash[s_node][t_node])
                                        dist_score.append(metabolite_score)
                                        accounted_metabolites[prod_reaction][t_node] = 1
                                    continue
                                else:
                                    # prod_reaction
                                    accounted[prod_reaction] = 1

                                reaction_score = 0 #score of each flux
                                if prod_reaction in reaction_folds:
                                    reaction_score = 1.0/(float(len_path)*denominator)*reaction_folds[prod_reaction]/len(path_hash[s_node][t_node])
                                    print prod_reaction, reaction_score, denominator, (float(len_path)*denominator)
                                dist_score.append(reaction_score)

                                #adding metabolite fold changes information, if available
                                metabolite_score = 0
                                if metabolite_folds and t_node in metabolite_folds:
                                    metabolite_score = 1.0/(float(len_path)*denominator)*metabolite_folds[t_node]/len(path_hash[s_node][t_node])
                                elif average_metabolite:
                                    metabolite_score = 1.0/(float(len_path)*denominator)*average_metabolite/len(path_hash[s_node][t_node])
                                dist_score.append(metabolite_score)
                                if metabolite_score != 0:
                                    accounted_metabolites[prod_reaction][t_node] = 1


                                if target + 1 != len(path):
                                    minus_reactions = set() #minusing consuming/producing reactions depending whether backward or forward
                                    node_consumers = B.successors(t_node) # consumption of target node, notice graph is reversed
                                    minus_reactions = set(node_consumers) - reaction_edges # minus reaction which are consumption of source not but are not used in shortest paths
                                    minus_reaction_score = 0
                                    #print "minus reactions", minus_reactions

                                    if minus_reactions and distance_global > 1: #minus consumption
                                        node_consumption_degree = float(len(minus_reactions))
                                        denominator_minus = denominator*node_consumption_degree
                                        for reaction in minus_reactions:
                                            if reaction in reaction_folds:
                                                minus_reaction_score = 1.0/(float(len_path)*denominator_minus)*reaction_folds[reaction]
                                                #print reaction, minus_reaction_score
                                                dist_score.append(-minus_reaction_score)

                                            metabolite_score = 0
                                            if metabolite_folds and t_node in metabolite_folds:
                                                metabolite_score = 1.0/(float(len_path)*denominator_minus)*metabolite_folds[t_node]/len(path_hash[s_node][t_node])
                                            elif average_metabolite:
                                                metabolite_score = 1.0/(float(len_path)*denominator_minus)*average_metabolite/len(path_hash[s_node][t_node])
                                            dist_score.append(-metabolite_score)

                            source = target
                            target += 1
            #consumption_score = 0.0
            #if cons1 and metabolite in consumption: #simply substracts consumption as in previous CoCCoA model
            #    consumption_score = consumption[metabolite]
            #scores[distance_global][metabolite] = sum(dist_score) - consumption_score
            scores[distance_global][metabolite] = sum(dist_score)
    return scores


def metabolites_genes(U, metabolites, consumption_genes, production_genes, degree=None):
    """For a given list of metabolites and degree returns genes (backward) for CoCCoA scores. By default calculates maximum possible degree.
    U: (DiGraph) graph of metabolites.
    metabolites: (list) list of metabolites.
    consumption/production: (dict) genes for production and consumption
    degree: (int) CoCCoA degree
    """
    if not isinstance(metabolites, list):
        raise TypeError("metabolites must be list")
    if not U.is_directed():
        raise TypeError("Graph must be directed")

    effective_degree = 0
    if degree >= 1:
        effective_degree = degree - 1 #effective degree since it is actually distance cuttoff
    if degree is None:
        effective_degree = None

    found_metabolites = list(set(U.nodes()) & set(metabolites))
    not_found = list(set(metabolites) - set(U.nodes()))

    if not_found:
        for m in not_found:
            print "No such metabolite {0}".format(m)

    if len(found_metabolites) == 0:
        raise TypeError("No metabolites found in the network")

    genes_results = collections.defaultdict(dict) #final results
    Urev = U.reverse() #reversing graph to get correct distances from source

    for metabolite in found_metabolites:
        path = nx.single_source_dijkstra_path_length(Urev, metabolite,
                                                     cutoff=effective_degree) #returns nodes and shortest path lengths until nodes, cutoff is maximal distance for search
        distances = {}
        for k, v in path.iteritems(): #reversing path lengths to get nodes which are at the distance from source metabolite distance[degree][metabolite]
            if v == 0 and degree == 0: #for zero degree
                distances[0] = distances.get(0, [])
                distances[0].append(k)
                break
            if v == 0: #for all the rest degrees
                distances[0] = distances.get(0, [])
                distances[0].append(k)
                distances[1] = distances.get(1, [])
                distances[1].append(k)
            elif v >= 1:
                n = v + 1 #
                distances[n] = distances.get(n, [])
                distances[n].append(k)
        score_genes = []
        for distance, mets in distances.iteritems():
            dist_genes = []
            if distance == 0: #special case
                cons_m = []
                if metabolite in consumption_genes:
                    cons_m = consumption_genes[metabolite]
                if cons_m:
                    dist_genes.extend(cons_m)
            elif distance == 1: #special case
                prod_m = []
                if metabolite in production_genes:
                    prod_m = production_genes[metabolite]
                if prod_m:
                    dist_genes.extend(prod_m)
            else:
                for m in mets:
                    cons_m = []
                    prod_m = []
                    if m in consumption_genes:
                        cons_m = consumption_genes[m]
                    if m in production_genes:
                        prod_m = production_genes[m]
                    if cons_m:
                        dist_genes.extend(cons_m)
                    if prod_m:
                        dist_genes.extend(prod_m)
            if dist_genes:
                score_genes.extend(dist_genes)
        genes_results[max(distances.keys())][metabolite] = list(set(score_genes))
    return genes_results


def metabolites_genes_forward(U, metabolites, consumption_genes, production_genes, degree=None):
    """For a given list of metabolites and degree returns genes for CoCCoA scores. By default calculates maximum possible degree.
    U: (DiGraph) graph of metabolites.
    metabolites: (list) list of metabolites.
    consumption/production: (dict) genes for production and consumption
    degree: (int) CoCCoA degree
    """
    if not isinstance(metabolites, list):
        raise TypeError("metabolites must be list")
    if not U.is_directed():
        raise TypeError("Graph must be directed")


    found_metabolites = list(set(U.nodes()) & set(metabolites))
    not_found = list(set(metabolites) - set(U.nodes()))

    if not_found:
        for m in not_found:
            print "No such metabolite {0}".format(m)

    if len(found_metabolites) == 0:
        raise TypeError("No metabolites found in the network")

    genes_results = collections.defaultdict(dict) #final results
    #Urev = U.reverse() #reversing graph to get correct distances from source

    for metabolite in found_metabolites:
        path_lengths = nx.single_source_dijkstra_path_length(U, metabolite,
                                                     cutoff=degree) #returns nodes and shortest path lengths until nodes, cutoff is maximal distance for search
        distances = {}
        for k, v in path_lengths.iteritems(): #reversing path lengths to get nodes which are at the distance from source metabolite distance[distance][metabolite]
            if v in distances:
                distances[v].append(k)
            else:
                distances[v] = [k]
        distances = collections.OrderedDict(sorted(distances.items()))
        score_genes = []
        for distance, mets in distances.iteritems():
            dist_genes = []
            if distance == 0:
                continue
            for m in mets:
                prod_m = []
                if m in production_genes:
                    prod_m = production_genes[m]
                if prod_m:
                    dist_genes.extend(prod_m)
            if dist_genes:
                score_genes.extend(dist_genes)
        genes_results[max(distances.keys())][metabolite] = list(set(score_genes))
    return genes_results


def metabolites_genes_backward(U, metabolites, consumption_genes, production_genes, degree=None):
    """For a given list of metabolites and degree returns genes for CoCCoA scores. By default calculates maximum possible degree.
    U: (DiGraph) graph of metabolites.
    metabolites: (list) list of metabolites.
    consumption/production: (dict) genes for production and consumption
    degree: (int) CoCCoA degree
    """
    if not isinstance(metabolites, list):
        raise TypeError("metabolites must be list")
    if not U.is_directed():
        raise TypeError("Graph must be directed")


    found_metabolites = list(set(U.nodes()) & set(metabolites))
    not_found = list(set(metabolites) - set(U.nodes()))

    if not_found:
        for m in not_found:
            print "No such metabolite {0}".format(m)

    if len(found_metabolites) == 0:
        raise TypeError("No metabolites found in the network")

    genes_results = collections.defaultdict(dict) #final results
    Urev = U.reverse() #reversing graph to get correct distances from source

    for metabolite in found_metabolites:
        path_lengths = nx.single_source_dijkstra_path_length(Urev, metabolite,
                                                     cutoff=degree) #returns nodes and shortest path lengths until nodes, cutoff is maximal distance for search
        distances = {}
        for k, v in path_lengths.iteritems(): #reversing path lengths to get nodes which are at the distance from source metabolite distance[distance][metabolite]
            if v in distances:
                distances[v].append(k)
            else:
                distances[v] = [k]
        distances = collections.OrderedDict(sorted(distances.items()))
        score_genes = []
        for distance, mets in distances.iteritems():
            dist_genes = []
            if distance == 0:
                continue
            for m in mets:
                cons_m = []
                if m in consumption_genes:
                    cons_m = consumption_genes[m]
                if cons_m:
                    dist_genes.extend(cons_m)
            if dist_genes:
                score_genes.extend(dist_genes)
        genes_results[max(distances.keys())][metabolite] = list(set(score_genes))
    return genes_results


def get_all_degrees_scores(U, metabolites, consumption, production):
    """Calculates all possible degree scores for given metabolite list. Returns dictionary of dictionaries
    results[degree][metabolite] = score.
    """
    scores = metabolites_scoring(U, metabolites, consumption, production, None) #determining maximum degree
    results = collections.defaultdict(dict)

    for degree, mets in scores.iteritems(): #going through all metabolites decreasing their degrees
        for m in mets:
            results[degree][m] = scores[degree][m]

        d = degree - 1
        while d >= 0:
            tmp_scores = metabolites_scoring(U, mets.keys(), consumption, production, d)
            for dd, s in tmp_scores.iteritems():
                for i in s:
                    results[dd][i] = tmp_scores[dd][i]
            d -= 1
    return results


def get_all_degrees_genes_backward(U, metabolites, consumption_genes, production_genes):
    """Gets all possible degree genes for given metabolite list. Returns dictionary of dictionaries
    results[degree][metabolite] = [genes].
    """
    genes = metabolites_genes_backward(U=U, metabolites=metabolites, consumption_genes=consumption_genes, production_genes=production_genes, degree=None) #determining maximum degree
    results = collections.defaultdict(dict)

    for degree, mets in genes.iteritems(): #going through all metabolites decreasing their degrees
        for m in mets:
            results[degree][m] = genes[degree][m]

        d = degree - 1
        while d >= 0:
            tmp_genes = metabolites_genes_backward(U=U, metabolites=mets.keys(), consumption_genes=consumption_genes, production_genes=production_genes, degree=d)
            for dd, s in tmp_genes.iteritems():
                for i in s:
                    results[dd][i] = tmp_genes[dd][i]
            d -= 1
    return results


def get_all_degrees_genes(U, metabolites, consumption_genes, production_genes):
    """Gets all possible degree genes for given metabolite list. Returns dictionary of dictionaries
    results[degree][metabolite] = [genes].
    """
    genes = metabolites_genes(U=U, metabolites=metabolites, consumption_genes=consumption_genes, production_genes=production_genes, degree=None) #determining maximum degree
    results = collections.defaultdict(dict)

    for degree, mets in genes.iteritems(): #going through all metabolites decreasing their degrees
        for m in mets:
            results[degree][m] = genes[degree][m]
        d = degree - 1
        while d >= 0:
            tmp_genes = metabolites_genes(U=U, metabolites=mets.keys(), consumption_genes=consumption_genes, production_genes=production_genes, degree=d)
            for dd, s in tmp_genes.iteritems():
                for i in s:
                    results[dd][i] = tmp_genes[dd][i]
            d -= 1
    return results


def get_all_degrees_genes_forward(U, metabolites, consumption_genes, production_genes):
    """Gets all possible degree genes for given metabolite list (forward direction). Returns dictionary of dictionaries
    results[degree][metabolite] = [genes].
    """
    genes = metabolites_genes_forward(U=U, metabolites=metabolites, consumption_genes=consumption_genes, production_genes=production_genes, degree=None) #determining maximum degree
    results = collections.defaultdict(dict)

    for degree, mets in genes.iteritems(): #going through all metabolites decreasing their degrees

        for m in mets:
            results[degree][m] = genes[degree][m]
        d = degree - 1
        while d >= 0:
            tmp_genes = metabolites_genes_forward(U=U, metabolites=mets.keys(), consumption_genes=consumption_genes, production_genes=production_genes, degree=d)
            for dd, s in tmp_genes.iteritems():
                for i in s:
                    results[dd][i] = tmp_genes[dd][i]
            d -= 1
    return results


def shuffle_expression(expression):
    genes = expression.keys()
    random.shuffle(genes)
    values = expression.values()
    shuffled = dict(zip(genes, values))
    return shuffled