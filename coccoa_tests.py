from unittest import TestCase
from coccoa import *
from routines import create_DiGraph, set_bounds, get_consumption_scores, get_production_scores, get_consumption_genes, \
    get_production_genes, metabolites_scoring, metabolite_scoring_forward, metabolite_scoring_bf, \
    metabolite_scoring_backward, get_all_degrees_genes_backward, \
    get_all_degrees_genes, get_all_degrees_genes_forward, shuffle_expression


class Tests(TestCase):
    def test_create_DiGraph_graph_vs_model(self):
        """Tests whether created graph is equal to model in terms of metabolites and reactions
        """

        model_text = """
-REACTIONS
R1: A + B -> 3 C
R2: B + C <-> 1 E
B3: B + C <-> 1 E
-CONSTRAINTS
R1[-100, 100]
R2[-100, 100]
-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        parser = BiooptParser()
        model = parser.parse(model_text)
        set_bounds(model, reverse=True)
        B = create_DiGraph(model, remove_nonfeas=False, remove_rev=False)
        self.assertTrue(bipartite.is_bipartite(B))

        met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
        react_nodes = set(B) - met_nodes


        ##checking metabolites
        metabolites_model = model.find_metabolites()
        met_names_model = [met.name for met in metabolites_model]
        self.assertEquals(met_nodes, set(met_names_model))

        ##checking reactions
        react_model = set([r.name for r in model.reactions])
        self.assertEquals(react_nodes, react_model)


    def test_create_DiGraph_bad_model1(self):
        model_text = """
-REACTIONS
R1: A + B <-> 3 C
R2: B + C <-> 1 E
R3: B + C <-> 1 E
R4: B + C <-> 1 E
-CONSTRAINTS
R1[-100, 100]
R2[-100, 100]

-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""

        parser = BiooptParser()
        model = parser.parse(model_text)
        set_bounds(model, reverse=True)
        self.assertRaises(RuntimeError, create_DiGraph, model=model, remove_rev=True, remove_nonfeas=True)

    def test_create_DiGraph_bad_model2(self):
        model_text = """
-REACTIONS
R1: A + B <-> 3 C
R2: B + C <-> 1 E
R3: B + C <-> 1 E
R4: B + C <-> 1 E
-CONSTRAINTS
R1[0, 0]
R2[0, 0]
R3[0, 0]
R4[0, 0]
-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
    """

        parser = BiooptParser()
        model = parser.parse(model_text)
        set_bounds(model, reverse=True)
        self.assertRaises(RuntimeError, create_DiGraph, model=model, remove_rev=True, remove_nonfeas=True)

    def test_create_DiGraph_bad_model3(self):
        model_text = """
-REACTIONS
R1: A + B <-> 3 C
A: B + C <-> 1 E
R3: B + C <-> 1 E
R4: B + C <-> 1 E
-CONSTRAINTS
R1[0, 100]
A[0, 100]

-EXTERNAL METABOLITES
E
-OBJ
A 1 1
-DESIGNOBJ
A R1 1
"""

        parser = BiooptParser()
        model = parser.parse(model_text)
        set_bounds(model, reverse=True)
        self.assertRaises(RuntimeError, create_DiGraph, model=model, remove_rev=True, remove_nonfeas=True)

    def test_create_DiGraph_good_model(self):
        model_text = """
-REACTIONS
R1: A + B <-> 3 C
R2: B + C <-> 1 E
R3: B + C <-> 1 E
R4: B + C <-> 1 E
-CONSTRAINTS
R1[-10, 0]
R2[0, 0]
R3[0, 10]
R4[0, 0]
-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
    """

        parser = BiooptParser()
        model = parser.parse(model_text)
        set_bounds(model, reverse=True)
        B = create_DiGraph(model=model, remove_rev=True, remove_nonfeas=True)

        G = nx.DiGraph()
        G.add_nodes_from(['A','B','C','E'], bipartite=0)
        G.add_nodes_from(['R1','R3'], bipartite=1)
        G.add_edges_from([('R1', 'A'), ('R1', 'B'), ('R3', 'E'), ('C', 'R1'), ('C', 'R3'), ('B', 'R3')])


        self.assertEquals(B.nodes(data=True), G.nodes(data=True))
        self.assertEquals(B.edges(), G.edges())


        #met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
        #react_nodes = set(B) - met_nodes


    def test_parse_minmax_bad_format(self):

        minmax_text = """first_line
#comment
R1\t10\t20
R2\t-10\t0
R3\t-10\t0
R4\t-10\t0\t2
"""
        self.assertRaises(RuntimeError, parse_minmax, minmax_text) # first line error
        self.assertRaises(RuntimeError, parse_minmax, minmax_text, skip_lines=1) # last line error

    def test_parse_minmax_tolerance(self):

        minmax_text = """firstline
R1\t-0.001\t20
R4\t-10\t1e-10
"""
        minmax_bounds = parse_minmax(minmax_text, skip_lines=1, tolerance=1e-7)
        self.assertEquals(minmax_bounds["R4"].ub, 0) #cheking tolerance towards 0
        self.assertEquals(minmax_bounds["R1"].lb, -0.001) #cheking tolerance towards 0

    def test_parse_minmax_duplicated_reactions(self):
        minmax_text = """first_line
#comment
R1\t10\t20
R2\t-10\t0
R2\t-10\t0
"""
        self.assertRaises(RuntimeError, parse_minmax, minmax_text, skip_lines=1, strict=True) # duplicate reactions strict

    def test_set_bounds_minmax(self):

        model_text = """
-REACTIONS
R1: A + B -> 3 C
R2: B + C <-> 1 E
B3: B + C <-> 1 E
-CONSTRAINTS
R1[-100, 100]
R2[-100, 100]
-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        parser = BiooptParser()
        model = parser.parse(model_text)
        minmax_bounds = []
        self.assertRaises(TypeError, set_bounds, model=model, minmax=minmax_bounds, reverse=True)   #min_max type check

        model = parser.parse(model_text)
        minmax_bounds = {}
        self.assertRaises(RuntimeError, set_bounds, model=model, minmax=minmax_bounds, reverse=True) #min_max empty

        model = parser.parse(model_text)
        minmax_bounds["R1"] = [1,2]
        self.assertRaises(TypeError, set_bounds, model=model, minmax=minmax_bounds, reverse=True)   #min_max without Bounds


    def test_set_bounds_reverse(self):
        model_text = """
-REACTIONS
R1: A + B -> 3 C
R2: B + C <-> 1 E
B3: B + C <-> 1 E
-CONSTRAINTS
R1[-100, 100]
R2[-100, 100]
-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        minmax_text = """
R1\t10\t20
R2\t-10\t-1
"""
        parser = BiooptParser()
        model = parser.parse(model_text)
        minmax_bounds = parse_minmax(minmax_text)
        react_before = model.find_reactions(minmax_bounds.keys())
        model_changed = set_bounds(copy.deepcopy(model), minmax_bounds, reverse=True)

        react_after = model_changed.find_reactions(minmax_bounds.keys())

        fwd = Direction.forward()
        rev = Direction.reversible()

        m_A = ReactionMember(Metabolite('A'), 1)
        m_B = ReactionMember(Metabolite('B'), 1)
        m_C = ReactionMember(Metabolite('C'), 1)
        m3_C = ReactionMember(Metabolite('C'), 3)
        m_E = ReactionMember(Metabolite('E', boundary=True), 1)

        R1 = Reaction('R1', m_A + m_B, m3_C, direction=fwd, bounds=Bounds(10, 20))
        R2 = Reaction('R2', m_E, m_B + m_C, direction=rev, bounds=Bounds(1, 10))

        self.assertEquals(react_after[1], R2) # checking if after reverse we have expected model with min_max provided

        model_text = """
-REACTIONS
R1: A + B -> 3 C
R2: B + C <-> 1 E
B3: B + C <-> 1 E
-CONSTRAINTS
R1[-100, 100]
R2[-10, -0]
-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""

        parser = BiooptParser()
        model = parser.parse(model_text)
        set_bounds(model, reverse=True)
        react_after = model.reactions

        m_B = ReactionMember(Metabolite('B'), 1)
        m_C = ReactionMember(Metabolite('C'), 1)
        m_E = ReactionMember(Metabolite('E', boundary=True), 1)
        R2 = Reaction('R2', m_E, m_B + m_C, direction=rev, bounds=Bounds(0, 10))
        self.assertEquals(react_after[1], R2) # checking if we have expected model after reversing reactions without minmax

    def test_metabolites_scoring(self):
        model_text = """
-REACTIONS
R1|R2: A -> C
R3: B -> C
R4: D <-> C
R5: C -> G
R6: G -> E
R7: G <-> F
R8: G -> F bla
-CONSTRAINTS
R1[0, 100]
R2[0, 100]
R3[0, 100]

-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        reaction2gene_text = """
R1|R2\tgene1\tgene2\tgene3
R3\tgene4
R4\tgene5\tgene6
R5\tgene7
R6\tgene8
R7\tgene9
"""
        expression_text = """
#comment
gene1\t5\t0.01
gene2\t5\t0.01
gene3\t5\t0.01
gene4\t-1\t0.01
gene5\t2\t0.01
gene6\t3\t0.01
gene7\t1\t0.01
gene8\t2\t0.01
gene9\t3\t0.01
"""
        minmax_text = """
#comment
R1\t0\t20
R4\t-10\t0
R2\t0\t10
R7\t0\t20
"""


        parser = BiooptParser()
        model = parser.parse(model_text)
        minmax_bounds = parse_minmax(minmax_text)
        set_bounds(model=model, minmax=minmax_bounds, reverse=True)
        B = create_DiGraph(model=model, remove_rev=True, remove_nonfeas=True)

        expression = read_expression(expression_text,sep="\t")
        reaction2gene = read_reaction2gene(reaction2gene_text)
        reaction_folds, reaction_genes = assign_expression(reaction2gene, expression)
        consumption = get_consumption_scores(B, reaction_folds)
        production = get_production_scores(B, reaction_folds)

        consumption_genes = get_consumption_genes(B, reaction_genes)
        production_genes = get_production_genes(B, reaction_genes)


        met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
        #print B.edges()
        U = nx.projected_graph(B, met_nodes) # projecting bipartite graph to make unipartite graph of metabolites
        Urev = U.reverse() #reversing graph to get correct distances from source
        #print U.nodes()
        #path = nx.single_source_dijkstra_path_length(Urev, "A", cutoff=4)
        #print path
        #cc = nx.connected_component_subgraphs(U.to_undirected())


        def get_maxdiameter(U):
            diameters = []
            H = U.to_undirected()
            subgraph_list = nx.connected_component_subgraphs(H)
            for G in subgraph_list:
                if len(G.nodes()) >= 1:
                    diameters.append(nx.diameter(G))
            if diameters:
                return max(diameters)
            else:
                return 0

        #print get_maxdiameter(U)

        metabolite = "F" #checking individual
        score_F_0 = metabolites_scoring(U=U, metabolites=[metabolite], consumption=consumption, production=production, degree=0)
        score_F_1 = metabolites_scoring(U=U, metabolites=[metabolite], consumption=consumption, production=production, degree=1)
        score_F_2 = metabolites_scoring(U=U, metabolites=[metabolite], consumption=consumption, production=production, degree=2)
        score_F_3 = metabolites_scoring(U=U, metabolites=[metabolite], consumption=consumption, production=production, degree=3)
        score_F_4 = metabolites_scoring(U=U, metabolites=[metabolite], consumption=consumption, production=production, degree=4)


        self.assertEquals(score_F_0[0]["F"], 0)
        self.assertEquals(score_F_1[1]["F"], -consumption["F"] + production["F"])
        self.assertEquals(score_F_2[2]["F"], -consumption["F"] + production["F"] + -consumption["G"] + production["G"])
        self.assertEquals(score_F_3[3]["F"], -consumption["F"] + production["F"] + -consumption["G"] + production["G"] + -consumption["C"] + production["C"])
        self.assertEquals(score_F_4[4]["F"], -consumption["F"] + production["F"] + -consumption["G"] + production["G"] + -consumption["C"] + production["C"] +
                                             0.5*(-consumption["A"] + production["A"] + -consumption["B"] + production["B"]) )

        #individual from the whole list
        metabolites = list(met_nodes)
        scores_0 = metabolites_scoring(U, metabolites, consumption, production, 0)
        scores_1 = metabolites_scoring(U, metabolites, consumption, production, 1)
        scores_2 = metabolites_scoring(U, metabolites, consumption, production, 2)
        scores_3 = metabolites_scoring(U, metabolites, consumption, production, 3)
        scores_4 = metabolites_scoring(U, metabolites, consumption, production, 4)

        self.assertEquals(scores_0[0]["C"], -1.5)
        self.assertEquals(scores_0[0]["F"], 0)
        self.assertEquals(scores_0[0]["E"], 0)
        self.assertEquals(scores_0[0]["G"], -2.5)
        self.assertEquals(scores_0[0]["D"], 0)
        self.assertEquals(scores_0[0]["A"], -5)
        self.assertEquals(scores_0[0]["B"], 1)

        self.assertEquals(scores_1[1]["A"], -5)
        self.assertEquals(scores_1[1]["B"], 1)
        self.assertEquals(scores_1[1]["C"], 0.5)
        self.assertEquals(scores_1[1]["D"], 2)
        self.assertEquals(scores_1[1]["E"], 2)
        self.assertEquals(scores_1[1]["F"], 3)
        self.assertEquals(scores_1[1]["G"], -1.5)

        self.assertNotIn("A", scores_2[2])
        self.assertNotIn("B", scores_2[2])
        self.assertEquals(scores_2[2]["C"], -1.5 + 2 -2)
        self.assertEquals(scores_2[2]["D"], 2 -1.5 + 2)
        self.assertEquals(scores_2[2]["E"], 2 -2.5 + 1)
        self.assertEquals(scores_2[2]["F"], 3 -2.5 + 1)
        self.assertEquals(scores_2[2]["G"], -2.5 + 1 -1.5 + 2)


        self.assertNotIn("A", scores_3[3])
        self.assertNotIn("B", scores_3[3])
        self.assertNotIn("C", scores_3[3])

        #print "cons", consumption
        #print "prod", production

        self.assertEquals(scores_3[3]["D"], -consumption["D"] + production["D"] + -consumption["C"] + production["C"] + 0.5*(-consumption["A"] + production["A"] + -consumption["B"] + production["B"]) )
        self.assertEquals(scores_3[3]["E"], -consumption["E"] + production["E"] + -consumption["G"] + production["G"] + -consumption["C"] + production["C"])
        self.assertEquals(scores_3[3]["F"], -consumption["F"] + production["F"] + -consumption["G"] + production["G"] + -consumption["C"] + production["C"])
        self.assertEquals(scores_3[3]["G"], -consumption["G"] + production["G"] + -consumption["C"] + production["C"] + 0.5*(-consumption["A"] + production["A"] + -consumption["B"] + production["B"]) )

        self.assertNotIn("A", scores_4[4])
        self.assertNotIn("B", scores_4[4])
        self.assertNotIn("C", scores_4[4])
        self.assertNotIn("D", scores_4[4])
        self.assertNotIn("G", scores_4[4])

        self.assertEquals(scores_4[4]["F"], -consumption["F"] + production["F"] + -consumption["G"] + production["G"] + -consumption["C"] + production["C"] +
                                             0.5*(-consumption["A"] + production["A"] + -consumption["B"] + production["B"]) )
        self.assertEquals(scores_4[4]["E"], -consumption["E"] + production["E"] + -consumption["G"] + production["G"] + -consumption["C"] + production["C"] +
                                             0.5*(-consumption["A"] + production["A"] + -consumption["B"] + production["B"]) )

        metabolite = "G" #checking individual if only 0 or 1st degree
        score_G_0 = metabolites_scoring(U, [metabolite], consumption, production, 0)
        score_G_1 = metabolites_scoring(U, [metabolite], consumption, production, 1)

        self.assertEquals(score_G_0[0]["G"], -consumption["G"])
        self.assertEquals(score_G_1[1]["G"], -consumption["G"] + production["G"])



        # self.assertIn("F", genes_F_1[1])
        #

        write_network_gml(B, "test.graphml", reaction_folds, reaction_genes, scoring=scores_2, degree=2)
        # metabolite = "C"
        # tmp = B.successors("C")
        # tmp.extend(B.predecessors("C"))
        # tmp.extend(metabolite)
        #
        # H = B.subgraph(tmp)
        #  #print reaction_folds
        #
        # for node in H.nodes():
        #      if node in reaction_genes:
        #          H.node[node]['genes'] = ';'.join(reaction_genes[node])
        # for edge in H.edges():
        #      if edge[0] in reaction_folds:
        #          H[edge[0]][edge[1]]['weight'] = reaction_folds[edge[0]]
        #      elif edge[1] in reaction_folds:
        #          H[edge[0]][edge[1]]['weight'] = reaction_folds[edge[1]]
        # print H.nodes(data=True)
        # print H.edges(data=True)
        #
        #nx.write_gml(H, "test.gml")

        #results = get_all_degrees_scores(U, list(met_nodes), consumption, production)
        #save_results(results, "test.txt")





    def test_read_expression_cutoff(self):
        expression_text="""
#comment
gene1\t10\t0.01
gene2\t5\t0.01
gene3\t-2\t0.2
"""
        expression = read_expression(expression_text, skip_lines=0, sep="\t", id=0, folds=1, pval=2,strict=False, cutoff=0.05)
        self.assertEquals(sorted(expression.keys()), ['gene1', 'gene2'])

    def test_read_expression_wrong_pvalue(self):
        expression_text = """
#comment
gene1\t10\t1.2
gene2\t5\t0.01
gene3\t-2\t0.2
"""
        self.assertRaises(ValueError, read_expression, text=expression_text, skip_lines=0, sep="\t", id=0, folds=1, pval=2, strict=False, cutoff=0.05)

    def test_read_expression_wrong_fold_change(self):
        expression_text="""
#comment
gene1\tbla\t1
gene2\t5\t0.01
gene3\t-2\t0.2
"""
        self.assertRaises(TypeError, read_expression, text=expression_text, skip_lines=0, sep="\t", id=0, folds=1, pval=2, strict=False, cutoff=0.05)

    def test_assign_expression(self):
        reaction2gene_text = """
R1|R2\tgene1\tgene2\tgene3
R3\tgene4
R4\tgene5\tgene6
"""
        expression_text = """
#comment
gene1\t10\t0.01
gene2\t5\t0.01
gene3\t-2\t0.01
gene4\t1\t0.01
gene5\t2\t0.01
gene6\t-3\t0.01
"""
        expression = read_expression(expression_text,sep="\t")
        reaction2gene = read_reaction2gene(reaction2gene_text)
        reaction_folds, r_gene_used = assign_expression(reaction2gene, expression)
        average = sum([expression['gene1'], expression['gene2'], expression['gene3']])/float(len([expression['gene1'], expression['gene2'], expression['gene3']]))
        minimum = min([expression['gene5'], expression['gene6']])
        self.assertEquals(reaction_folds['R1|R2'], average)
        self.assertEquals(reaction_folds['R3'], 1)
        self.assertEquals(reaction_folds['R4'], minimum)


    def test_shuffle_expression(self):
        expression_text = """
#comment
gene1\t10\t0.01
gene2\t5\t0.01
gene3\t-2\t0.01
gene4\t1\t0.01
gene5\t2\t0.01
gene6\t-3\t0.01
"""
        expression = read_expression(expression_text,sep="\t")
        a = shuffle_expression(expression)
        self.assertNotEqual(expression, a)



    def test_metabolite_scoring2(self):

        model_text = """
-REACTIONS
R1|R2: A -> C
R3: B -> C
R4: D <-> C
R5: C -> G
R6: G -> E
R7: G <-> F
-CONSTRAINTS
R1[0, 100]
R2[0, 100]
R3[0, 100]

-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        reaction2gene_text = """
R1|R2\tgene1\tgene2\tgene3
R3\tgene4
R4\tgene5\tgene6
R5\tgene7
R6\tgene8
R7\tgene9
"""
        expression_text = """
#comment
gene1\t5\t0.01
gene2\t5\t0.01
gene3\t5\t0.01
gene4\t-1\t0.01
gene5\t2\t0.01
gene6\t3\t0.01
gene7\t1\t0.01
gene8\t2\t0.01
gene9\t3\t0.01
"""
        minmax_text = """
#comment
R1\t0\t20
R4\t-10\t0
R2\t0\t10
R7\t0\t20
"""

        parser = BiooptParser()
        model = parser.parse(model_text)
        minmax_bounds = parse_minmax(minmax_text)
        set_bounds(model=model, minmax=minmax_bounds, reverse=True)
        B = create_DiGraph(model=model, remove_rev=True, remove_nonfeas=True)

        expression = read_expression(expression_text,sep="\t")
        reaction2gene = read_reaction2gene(reaction2gene_text)
        reaction_folds, reaction_genes = assign_expression(reaction2gene, expression)
        consumption = get_consumption_scores(B, reaction_folds)
        production = get_production_scores(B, reaction_folds)

        #consumption_genes = get_consumption_genes(B, reaction_genes)
        #production_genes = get_production_genes(B, reaction_genes)


        met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
        #print B.edges()
        U = nx.projected_graph(B, met_nodes) # projecting bipartite graph to make unipartite graph of metabolites
        Urev = U.reverse() #reversing graph to get correct distances from source
        #print U.nodes()
        path = nx.single_source_dijkstra_path_length(Urev, "C", cutoff=1)
        #print path




        # def get_consumption_scores(B, reaction_folds):
        #     """For each metabolite in B returns dictionary average reaction fold-changes from production side
        #     """
        #     scores = {}
        #     met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
        #     for metabolite in met_nodes:
        #         tmp_score = []
        #         for reaction in B.successors(metabolite):
        #             if reaction in reaction_folds:
        #                 tmp_score.append(reaction_folds[reaction])
        #         if tmp_score:
        #             scores[metabolite] = float(sum(tmp_score)) / float(len(tmp_score))
        #         else:
        #             scores[metabolite] = float(0)
        #     if not scores:
        #         raise RuntimeError("None of production scores were computed")
        #     return scores


        # def get_production_scores(B, consumption):
        #     """For each metabolite in B returns dictionary average reaction fold-changes from consumption side
        #     """
        #
        #     scores = {}
        #     met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
        #     U = nx.projected_graph(B, met_nodes)
        #     Urev = U.reverse()
        #
        #     for metabolite in met_nodes:
        #         path = nx.single_source_dijkstra_path_length(Urev, metabolite, cutoff=1)
        #         n_score = []
        #         for neighbour in path:
        #             if path[neighbour] == 1:
        #                 if neighbour in consumption and consumption[neighbour] != 0:
        #                     n_score.append(consumption[neighbour])
        #         if n_score:
        #             scores[metabolite] = float(sum(n_score)) / float(len(n_score))
        #         else:
        #             scores[metabolite] = float(0)
        #
        #     if not scores:
        #         raise RuntimeError("None of production scores were computed")
        #     return scores

        consumption = get_consumption_scores(B, reaction_folds)
        production = get_production_scores(B, reaction_folds)

        def metabolites_scoring_forward(U, metabolites, consumption, production, degree=None):
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
            #Urev = U.reverse() #reversing graph to get correct distances from source

            for metabolite in found_metabolites:
                path = nx.single_source_dijkstra_path_length(U, metabolite,
                                                          cutoff=effective_degree) #returns nodes and shortest path lengths until nodes, cutoff is maximal distance for search
                distances = {}
                print metabolite, path
                for k, v in path.iteritems(): #reversing path lengths to get nodes which are at the distance from source metabolite distance[degree][metabolite]
                    if v == 0 and degree == 0: #for zero degree
                        break
                    if v == 0: #for all the rest degrees
                        distances[1] = [k]
                    elif v >= 1:
                        n = v + 1 #
                        if n in distances:
                            distances[n].append(k)
                        else:
                            distances[n] = [k]

                print "forward path from function", path, distances
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


        def metabolites_scoring_backward(U, metabolites, consumption, production, degree=None):
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
                print metabolite, path
                for k, v in path.iteritems(): #reversing path lengths to get nodes which are at the distance from source metabolite distance[degree][metabolite]
                    if v == 0 and degree == 0: #for zero degree
                        break
                    if v == 0: #for all the rest degrees
                        continue
                    elif v >= 1:
                        n = v + 1 #
                        if n in distances:
                            distances[n].append(k)
                        else:
                            distances[n] = [k]

                print "backward path from function", path, distances

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

        # scores_f = metabolites_scoring_forward(U, ["C"], consumption=consumption, production=production, degree=1)
        # scores_b = metabolites_scoring_backward(U, ["C"], consumption=consumption, production=production, degree=2)
        #
        # path_back    = nx.single_source_dijkstra_path_length(Urev, "C", cutoff=2)
        # path_forward = nx.single_source_dijkstra_path_length(U, "C", cutoff=2)
        #
        # print "Backward:", path_back
        # print "Forward:", path_forward
        # print "Scores forward", scores_f
        # print "Scores_backward", scores_b

    def test_metabolites_scoring_backward(self):
        model_text = """
-REACTIONS
R1|R2: A -> C
R3: B -> C
R4: D <-> C
R5: C -> G
R6: G -> E
R7: G <-> F
R8: I + H -> G
R9: C -> H
R10: H -> Z

-CONSTRAINTS
R1[0, 100]
R2[0, 100]
R3[0, 100]

-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        reaction2gene_text = """
R1|R2\tgene1\tgene2\tgene3
R3\tgene4
R4\tgene5\tgene6
R5\tgene7
R6\tgene8
R7\tgene9
R8\tgene10
R9\tgene11
R10\tgene12
R11\tgene13
R12\tgene14
"""
        expression_text = """
#comment
gene1\t5\t0.01
gene2\t5\t0.01
gene3\t5\t0.01
gene4\t-1\t0.01
gene5\t1\t0.01
gene6\t3\t0.01
gene7\t3\t0.01
gene8\t2\t0.01
gene9\t0.5\t0.01
gene10\t-2\t0.01
gene11\t5\t0.01
gene12\t2\t0.01
gene13\t-0.5\t0.01
gene14\t1\t0.01

"""
        minmax_text = """
#comment
R1\t0\t20
R4\t-10\t0
R2\t0\t10
R7\t0\t20
"""

        parser = BiooptParser()
        model = parser.parse(model_text)
        minmax_bounds = parse_minmax(minmax_text)
        set_bounds(model=model, minmax=minmax_bounds, reverse=True)
        B = create_DiGraph(model=model, remove_rev=True, remove_nonfeas=True)

        expression = read_expression(expression_text,sep="\t")
        reaction2gene = read_reaction2gene(reaction2gene_text)
        reaction_folds, reaction_genes = assign_expression(reaction2gene, expression)

        metabolite = "G"
        scores = metabolite_scoring_backward(B=B, metabolites=[metabolite], reaction_folds=reaction_folds, degree=2)

        self.assertEquals(scores[0][metabolite], 0)
        self.assertEquals(scores[1][metabolite], 1.0/2.0*reaction_folds["R5"] + 1.0/2.0*reaction_folds["R8"])
        #self.assertEquals(scores[2][metabolite], 1.0/4.0*reaction_folds["R5"] + 1.0/2.0*reaction_folds["R8"] + 1.0/8.0*reaction_folds["R1|R2"] + 1.0/8.0*reaction_folds["R3"] -1.0/4.0*reaction_folds["R4"] )
        self.assertEquals(scores[2][metabolite], 1.0/4.0*reaction_folds["R5"] + 1.0/2.0*reaction_folds["R8"] + 1.0/8.0*reaction_folds["R1|R2"] + 1.0/8.0*reaction_folds["R3"] -1.0/8.0*reaction_folds["R4"] -1.0/8.0*reaction_folds["R9"] )




    def test_metabolites_scoring_forward(self):
        model_text = """
-REACTIONS
R1|R2: A -> C
R3: B -> C
R4: D <-> C
R5: C -> G
R6: G -> E
R7: G <-> F
R8: H -> G
R9: H -> I


-CONSTRAINTS
R1[0, 100]
R2[0, 100]
R3[0, 100]

-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        reaction2gene_text = """
R1|R2\tgene1\tgene2\tgene3
R3\tgene4
R4\tgene5\tgene6
R5\tgene7
R6\tgene8
R7\tgene9
R8\tgene10
R9\tgene11
R10\tgene12
R11\tgene13
R12\tgene14
"""
        expression_text = """
#comment
gene1\t5\t0.01
gene2\t5\t0.01
gene3\t5\t0.01
gene4\t-1\t0.01
gene5\t1\t0.01
gene6\t3\t0.01
gene7\t3\t0.01
gene8\t2\t0.01
gene9\t0.5\t0.01
gene10\t-2\t0.01
gene11\t1\t0.01
gene12\t2\t0.01
gene13\t-0.5\t0.01
gene14\t1\t0.01

"""
        minmax_text = """
#comment
R1\t0\t20
R4\t-10\t0
R2\t0\t10
R7\t0\t20
"""

        parser = BiooptParser()
        model = parser.parse(model_text)
        minmax_bounds = parse_minmax(minmax_text)
        set_bounds(model=model, minmax=minmax_bounds, reverse=True)
        B = create_DiGraph(model=model, remove_rev=True, remove_nonfeas=True)

        expression = read_expression(expression_text,sep="\t")
        reaction2gene = read_reaction2gene(reaction2gene_text)
        reaction_folds, reaction_genes = assign_expression(reaction2gene, expression)

        consumption = get_consumption_scores(B, reaction_folds)
        production = get_production_scores(B, reaction_folds)

        metabolite = "C"
        scores = metabolite_scoring_forward(B, [metabolite], reaction_folds, degree=2)

        self.assertEquals(scores[0][metabolite], 0)
        self.assertEquals(scores[1][metabolite], 1.0/2.0*reaction_folds["R4"] + 1.0/2.0*reaction_folds["R5"])
        #self.assertEquals(scores[2][metabolite], 1.0/2.0*reaction_folds["R4"] + 1.0/4.0*reaction_folds["R5"] + 1.0/8.0*reaction_folds["R6"] + 1.0/8.0*reaction_folds["R7"])
        self.assertEquals(scores[2][metabolite], 1.0/2.0*reaction_folds["R4"] + 1.0/4.0*reaction_folds["R5"] + 1.0/8.0*reaction_folds["R6"] + 1.0/8.0*reaction_folds["R7"] -1.0/4.0*reaction_folds["R8"] )



    def test_metabolites_scoring_backward_met(self):
        model_text = """
-REACTIONS
R1|R2: A -> C
R3: B -> C
R4: D <-> C
R5: C -> G
R6: G -> E
R7: G <-> F
R8: I + H -> G
R9: H -> Z

-CONSTRAINTS
R1[0, 100]
R2[0, 100]
R3[0, 100]

-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        reaction2gene_text = """
R1|R2\tgene1\tgene2\tgene3
R3\tgene4
R4\tgene5\tgene6
R5\tgene7
R6\tgene8
R7\tgene9
R8\tgene10
R9\tgene11
R10\tgene12
R11\tgene13
R12\tgene14
"""
        metabolites_text = """
#comment
A\t1\t0.01
B\t2\t0.01
C\t3\t0.01
D\t-2.5\t0.01
E\t-4\t0.01
F\t0.5\t0.01
G\t1\t0.01
"""
        minmax_text = """
#comment
R1\t0\t20
R4\t-10\t0
R2\t0\t10
R7\t0\t20
"""

        expression_text = """
#comment
gene1\t5\t0.01
gene2\t5\t0.01
gene3\t5\t0.01
gene4\t-1\t0.01
gene5\t1\t0.01
gene6\t3\t0.01
gene7\t3\t0.01
gene8\t2\t0.01
gene9\t0.5\t0.01
gene10\t-2\t0.01
gene11\t1\t0.01
gene12\t2\t0.01
gene13\t-0.5\t0.01
gene14\t1\t0.01
gene15\t1\t0.01
"""


        parser = BiooptParser()
        model = parser.parse(model_text)
        minmax_bounds = parse_minmax(minmax_text)
        set_bounds(model=model, minmax=minmax_bounds, reverse=True)
        B = create_DiGraph(model=model, remove_rev=True, remove_nonfeas=True)

        expression = read_expression(expression_text,sep="\t")
        metabolite_folds = read_metabolites(metabolites_text, sep="\t")

        reaction2gene = read_reaction2gene(reaction2gene_text)
        reaction_folds, reaction_genes = assign_expression(reaction2gene, expression)
        metabolite = "G"
        avg_metabolite = sum(metabolite_folds.values())/float(len(metabolite_folds.values()))
        scores = metabolite_scoring_backward(B, [metabolite], reaction_folds=reaction_folds, metabolite_folds=metabolite_folds, degree=2)

        self.assertEquals(round(scores[1][metabolite],5), round(1.0/2.0*reaction_folds["R5"] + 1.0/2.0*metabolite_folds["C"] + 1.0/2.0*reaction_folds["R8"] + 1.0/2.0*avg_metabolite + 1.0/2.0*avg_metabolite,5))
        self.assertEquals(scores[2][metabolite], -1/4.0*reaction_folds["R4"] - 1/4.0*metabolite_folds["C"] + 1.0/4.0*(reaction_folds["R5"] + metabolite_folds["C"]) + 1.0/2.0*(reaction_folds["R8"]+avg_metabolite + avg_metabolite) + 1.0/8.0*(reaction_folds["R1|R2"] + metabolite_folds["A"]) + 1.0/8.0*(reaction_folds["R3"] + metabolite_folds["B"]) )


    def test_metabolites_scoring_forward_met(self):
        model_text = """
-REACTIONS
R1|R2: A -> C
R3: B -> C
R4: D <-> C
R5: C -> G
R6: G -> E
R7: G <-> F
R8: H -> G



-CONSTRAINTS
R1[0, 100]
R2[0, 100]
R3[0, 100]

-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        reaction2gene_text = """
R1|R2\tgene1\tgene2\tgene3
R3\tgene4
R4\tgene5\tgene6
R5\tgene7
R6\tgene8
R7\tgene9
R8\tgene10
R9\tgene11
R10\tgene12
R11\tgene13
R12\tgene14
"""
        metabolites_text = """
#comment
A\t1\t0.01
B\t2\t0.01
C\t3\t0.01
D\t-2.5\t0.01
E\t-4\t0.01
F\t0.5\t0.01
G\t5\t0.01
"""
        minmax_text = """
#comment
R1\t0\t20
R4\t-10\t0
R2\t0\t10
R7\t0\t20
"""

        expression_text = """
#comment
gene1\t5\t0.01
gene2\t5\t0.01
gene3\t5\t0.01
gene4\t-1\t0.01
gene5\t1\t0.01
gene6\t3\t0.01
gene7\t3\t0.01
gene8\t2\t0.01
gene9\t0.5\t0.01
gene10\t-2\t0.01
gene11\t1\t0.01
gene12\t2\t0.01
gene13\t-0.5\t0.01
gene14\t1\t0.01
"""


        parser = BiooptParser()
        model = parser.parse(model_text)
        minmax_bounds = parse_minmax(minmax_text)
        set_bounds(model=model, minmax=minmax_bounds, reverse=True)
        B = create_DiGraph(model=model, remove_rev=True, remove_nonfeas=True)

        expression = read_expression(expression_text,sep="\t")
        metabolite_folds = read_metabolites(metabolites_text, sep="\t")

        reaction2gene = read_reaction2gene(reaction2gene_text)
        reaction_folds, reaction_genes = assign_expression(reaction2gene, expression)
        metabolite = "C"
        avg_metabolite = sum(metabolite_folds.values())/float(len(metabolite_folds.values()))
        scores = metabolite_scoring_forward(B, [metabolite], reaction_folds=reaction_folds, degree=2)
        self.assertEquals(scores[1][metabolite], 1/2.0*(reaction_folds["R5"] + reaction_folds["R4"]))
        self.assertEquals(scores[2][metabolite], 1/4.0*(reaction_folds["R5"] + reaction_folds["R4"]) - 1.0/4.0*reaction_folds["R8"] + 1/4.0*reaction_folds["R4"] + 1.0/8.0*(reaction_folds["R6"] + reaction_folds["R7"]) )

        scores = metabolite_scoring_forward(B, metabolites=[metabolite], reaction_folds=reaction_folds, metabolite_folds=metabolite_folds, degree=2)
        self.assertEquals(scores[1][metabolite], 1/2.0*(reaction_folds["R5"] + reaction_folds["R4"]))
        self.assertEquals(scores[2][metabolite], 1/4.0*(reaction_folds["R5"] + reaction_folds["R4"]) - 1.0/4.0*(reaction_folds["R8"]+avg_metabolite) + 1/4.0*reaction_folds["R4"] + 1.0/8.0*(reaction_folds["R6"] + reaction_folds["R7"] + 2*metabolite_folds["G"]) )

    def test_metabolite_genes(self):

        model_text = """
-REACTIONS
R1|R2: A -> C
R3: B -> C
R4: D <-> C
R5: C -> G
R6: G -> E
R7: G <-> F
-CONSTRAINTS
R1[0, 100]
R2[0, 100]
R3[0, 100]

-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        reaction2gene_text = """
R1|R2\tgene1\tgene2\tgene3
R3\tgene4
R4\tgene5\tgene6
R5\tgene7
R6\tgene8
R7\tgene9
"""
        expression_text = """
#comment
gene1\t5\t0.01
gene2\t5\t0.01
gene3\t5\t0.01
gene4\t-1\t0.01
gene5\t2\t0.01
gene6\t3\t0.01
gene7\t1\t0.01
gene8\t2\t0.01
gene9\t3\t0.01
"""
        minmax_text = """
#comment
R1\t0\t20
R4\t-10\t0
R2\t0\t10
R7\t0\t20
"""


        parser = BiooptParser()
        model = parser.parse(model_text)
        minmax_bounds = parse_minmax(minmax_text)
        set_bounds(model=model, minmax=minmax_bounds, reverse=True)
        B = create_DiGraph(model=model, remove_rev=True, remove_nonfeas=True)

        expression = read_expression(expression_text,sep="\t")
        reaction2gene = read_reaction2gene(reaction2gene_text)
        reaction_folds, reaction_genes = assign_expression(reaction2gene, expression)
        consumption = get_consumption_scores(B, reaction_folds)
        production = get_production_scores(B, reaction_folds)

        consumption_genes = get_consumption_genes(B, reaction_genes)
        production_genes = get_production_genes(B, reaction_genes)

        print "Production:"
        print production_genes

        print "Consumption::"
        print consumption_genes

        met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
        #print B.edges()
        U = nx.projected_graph(B, met_nodes) # projecting bipartite graph to make unipartite graph of metabolites

        metabolite = "G"
        results_genes = get_all_degrees_genes(U, [metabolite], consumption_genes=consumption_genes, production_genes=production_genes)

        self.assertEquals(sorted(results_genes[0]["G"]), ["gene8", "gene9"])
        self.assertEquals(sorted(results_genes[1]["G"]), ["gene7","gene8", "gene9"])
        self.assertEquals(sorted(results_genes[2]["G"]), ["gene1", "gene2", "gene3", "gene4", "gene5", "gene7", "gene8", "gene9"])


        results_genes_backward = get_all_degrees_genes_backward(U, [metabolite], consumption_genes=consumption_genes, production_genes=production_genes)
        self.assertEquals(sorted(results_genes_backward[0]["G"]), [])
        self.assertEquals(sorted(results_genes_backward[1]["G"]), ["gene5", "gene7"])
        self.assertEquals(sorted(results_genes_backward[2]["G"]), ["gene1", "gene2", "gene3", "gene4", "gene5", "gene7"])

        results_genes_forward = get_all_degrees_genes_forward(U, ["A"], consumption_genes=consumption_genes, production_genes=production_genes)
        self.assertEquals(results_genes_forward[0]["A"], [])
        self.assertEquals(sorted(results_genes_forward[1]["A"]), ["gene1","gene2","gene3", "gene4"])
        self.assertEquals(sorted(results_genes_forward[2]["A"]), ["gene1","gene2","gene3", "gene4","gene5","gene7"])


    def test_metabolites_scoring_backward_test(self):
        model_text = """
-REACTIONS
R1|R2: A -> C
R1: J -> A
R3: B -> C
R4: D <-> C
R5: C -> G
R6: G -> E
R7: G <-> F
R8: I + H -> G
R9: C -> H
R10: H -> Z
R11: C -> G
R12: K -> G

-CONSTRAINTS
R1[0, 100]
R2[0, 100]
R3[0, 100]

-EXTERNAL METABOLITES
E
-OBJ
R2 1 1
-DESIGNOBJ
R2 R1 1
"""
        reaction2gene_text = """
R1|R2\tgene1\tgene2\tgene3
R1\tgene15
R3\tgene4
R4\tgene5\tgene6
R5\tgene7
R6\tgene8
R7\tgene9
R8\tgene10
R9\tgene11
R10\tgene12
R11\tgene13
R12\tgene14
"""
        expression_text = """
#comment
gene1\t5\t0.01
gene2\t5\t0.01
gene3\t5\t0.01
gene4\t-1\t0.01
gene5\t1\t0.01
gene6\t3\t0.01
gene7\t3\t0.01
gene8\t2\t0.01
gene9\t0.5\t0.01
gene10\t-2\t0.01
gene11\t5\t0.01
gene12\t2\t0.01
gene13\t-0.5\t0.01
gene14\t1\t0.01

"""
        minmax_text = """
#comment
R1\t0\t20
R4\t-10\t0
R2\t0\t10
R7\t0\t20
"""

        metabolites_text = """
#comment
A\t1\t0.01
B\t2\t0.01
C\t3\t0.01
D\t-2.5\t0.01
E\t-4\t0.01
F\t0.5\t0.01
G\t5\t0.01
"""

        model_parser = BiooptParser()
        model = model_parser.parse(model_text)
        warnings.simplefilter("ignore")
        #infile_name = "Z:/Alex/transcriptome_metabolome/results/2013-07-07/tmp/bounds2_fixedAA.txt"
        #model = model_parser.parse_file(infile_name)
        #minmax_file = "Z:/Alex/transcriptome_metabolome/results/2013-07-07/tmp/minmax2AA.txt"
        minmax_bounds = parse_minmax(minmax_text)
        set_bounds(model=model, minmax=minmax_bounds, reverse=True)
        B = create_DiGraph(model=model, remove_rev=True, remove_nonfeas=True)
        met_nodes = set(n for n, d in B.nodes(data=True) if d['bipartite'] == 0)
        U = nx.projected_graph(B, met_nodes) # projecting bipartite graph to make unipartite graph of metabolites
        Urev = U.reverse()

        expression = read_expression(expression_text,sep="\t")
        metabolite_folds = read_metabolites(metabolites_text, sep="\t")

        reaction2gene = read_reaction2gene(reaction2gene_text)
        reaction_folds, reaction_genes = assign_expression(reaction2gene, expression)

        metabolite = "G"

        print "########TESTING#################\n\n\n"
        scores = metabolite_scoring_bf(B=B, metabolites=[metabolite], reaction_folds=reaction_folds, degree=2)

        self.assertEquals(scores[1][metabolite], 1.0/4.0*reaction_folds["R5"] + 1.0/4.0*reaction_folds["R11"] + 1.0/4.0*reaction_folds["R12"]
                                                 + 1/2.0*1.0/4.0*(reaction_folds["R8"]+reaction_folds["R9"]) - 1/8.0*reaction_folds["R10"])

        self.assertEquals(round(scores[2][metabolite],5), round(3/7.0*1/2.0*reaction_folds["R3"] + 3/7.0*1/2.0*reaction_folds["R1|R2"]
                                                 + 2/7.0*1/2.0*reaction_folds["R5"] + 2/7.0*1/2.0*reaction_folds["R11"]
                                                 + 2/7.0*1/3.0*reaction_folds["R8"] + 2/7.0*1/3.0*reaction_folds["R9"]
                                                 + 1/7.0*1/1.0*reaction_folds["R12"] - 2*2/7.0*1.0/2.0*reaction_folds["R4"]
                                                 - 2/7.0*1.0/3.0*reaction_folds["R4"]- 2/7.0*1.0/3.0*reaction_folds["R10"],5))

        scores_met = metabolite_scoring_bf(B=B, metabolites=[metabolite], reaction_folds=reaction_folds,  metabolite_folds=metabolite_folds, degree=2)
        avg_metabolite = sum(metabolite_folds.values())/float(len(metabolite_folds.values()))
        self.assertEquals(round(scores_met[1][metabolite],5), round(1/4.0*reaction_folds["R5"] + 1/4.0*reaction_folds["R11"] + 1/4.0*reaction_folds["R12"]
                                                  + 1/2.0*1/4.0*(reaction_folds["R9"]+reaction_folds["R8"]) - 1/8.0*reaction_folds["R10"]
                                                  + 1/4.0*metabolite_folds["C"] + 1/4.0*metabolite_folds["C"] + 1/4.0*avg_metabolite
                                                  + 1/2.0*1/4.0*metabolite_folds["C"] + 2*1/2.0*1/4.0*avg_metabolite*1/2.0 - 1/8.0*avg_metabolite, 5))

        self.assertEquals(round(scores_met[2][metabolite],5), round(3/7.0*1/2.0*reaction_folds["R3"] + 3/7.0*1/2.0*reaction_folds["R1|R2"]
                                                  + 2/7.0*1/2.0*reaction_folds["R5"] + 2/7.0*1/2.0*reaction_folds["R11"]
                                                  + 2/7.0*1/3.0*reaction_folds["R8"] + 2/7.0*1/3.0*reaction_folds["R9"]
                                                  + 1/7.0*1/1.0*reaction_folds["R12"]
                                                  - 2*2/7.0*1.0/2.0*reaction_folds["R4"] - 2/7.0*1.0/3.0*reaction_folds["R4"]
                                                  - 2/7.0*1.0/3.0*reaction_folds["R10"]
                                                  + 3/7.0*1/2.0*metabolite_folds["B"] + 3/7.0*1/2.0*metabolite_folds["A"]
                                                  + 2/7.0*1/2.0*metabolite_folds["C"] + 2/7.0*1/2.0*metabolite_folds["C"]
                                                  + 2.0*2/7.0*1/3.0*avg_metabolite/2.0 + 2/7.0*1/3.0*metabolite_folds["C"]
                                                  + 1/7.0*1/1.0*avg_metabolite
                                                  - 2*2/7.0*1.0/2.0*metabolite_folds["C"] - 2/7.0*1.0/3.0*metabolite_folds["C"]
                                                  - 2/7.0*1.0/3.0*avg_metabolite,5) )
        metabolite = "C"
        scores_forward = metabolite_scoring_bf(B=B, metabolites=[metabolite], reaction_folds=reaction_folds, degree=2, dir="f")

        self.assertEquals(scores_forward[1][metabolite], 1.0/4.0*reaction_folds["R4"] + 1.0/4.0*reaction_folds["R5"] + 1.0/4.0*reaction_folds["R11"]
                                                       + 1/2.0*1.0/4.0*(reaction_folds["R8"]+reaction_folds["R9"]))

        self.assertEquals(scores_forward[2][metabolite], 3.0/8.0*1/2.0*reaction_folds["R6"] + 3.0/8.0*1/2.0*reaction_folds["R7"]
                                                       + 3.0/8.0*1/2.0*reaction_folds["R9"] + 2.0/8.0*1/2.0*reaction_folds["R5"]
                                                       + 2.0/8.0*1/2.0*reaction_folds["R11"] + 2.0/8.0*1/3.0*reaction_folds["R8"]
                                                       + 1.0/8.0*1*reaction_folds["R4"] + 1.0/8.0*1/2.0*reaction_folds["R10"]
                                                       - 2.0/8.0*1/3.0*reaction_folds["R12"] - 2.0/8.0*1/2.0*reaction_folds["R12"]
                                                       - 2.0/8.0*1/2.0*reaction_folds["R12"])

        metabolite = "C"
        scores_forward_met = metabolite_scoring_bf(B=B, metabolites=[metabolite], reaction_folds=reaction_folds, metabolite_folds=metabolite_folds, degree=2, dir="f")
        self.assertEquals(round(scores_forward_met[1][metabolite], 5), round(1.0/4.0*reaction_folds["R4"] + 1.0/4.0*reaction_folds["R5"] + 1.0/4.0*reaction_folds["R11"]
                                                             + 1/2.0*1.0/4.0*(reaction_folds["R9"] + reaction_folds["R8"]) + 2*1/2.0*1.0/4.0*avg_metabolite/2.0, 5))

        self.assertEquals(scores_forward_met[2][metabolite],   3.0/8.0*1/2.0*reaction_folds["R6"] + 3.0/8.0*1/2.0*metabolite_folds["G"]
                                                         + 3.0/8.0*1/2.0*reaction_folds["R7"] + 3.0/8.0*1/2.0*metabolite_folds["G"]
                                                         + 3.0/8.0*1/2.0*reaction_folds["R9"]
                                                         + 2.0/8.0*1/2.0*reaction_folds["R5"]
                                                         + 2.0/8.0*1/2.0*reaction_folds["R11"]
                                                         + 2.0/8.0*1/3.0*reaction_folds["R8"] + 2*2.0/8.0*1/3.0*avg_metabolite/2.0
                                                         + 1.0/8.0*1*reaction_folds["R4"]
                                                         + 1.0/8.0*1/2.0*reaction_folds["R10"] + 1.0/8.0*1/2.0*avg_metabolite
                                                         - 2.0/8.0*1/3.0*reaction_folds["R12"] - 2.0/8.0*1/3.0*avg_metabolite
                                                         - 2.0/8.0*1/2.0*reaction_folds["R12"] - 2.0/8.0*1/2.0*avg_metabolite
                                                         - 2.0/8.0*1/2.0*reaction_folds["R12"] - 2.0/8.0*1/2.0*avg_metabolite)

