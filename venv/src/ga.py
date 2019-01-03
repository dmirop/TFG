#!/usr/bin/python

import pool
import chromosome as chroms
from random import choices, choice, sample
from math import sqrt
import heapq
import copy as cp
import time
import multiprocessing as mp


def create_file_name():
    """
    Creates a string to identify de file name
    :return: a string that contains the full date
    """
    readable = (str(time.ctime(time.time()))).replace(":", "_")
    return readable


class GeneticAlgorithm:
    """
    Defines a generic Genetic Algorithm
    """
    def __init__(self, pool_size=100, p_cross=0.8, p_muta=0.05, elitism=True, endogamy=True, max_gen=5000, max_change=1500):
        self._pool = pool.Pool()
        self._pool_size = pool_size
        self._p_cross = p_cross
        self._p_muta = p_muta
        self._elitism = elitism
        self._endogamy = endogamy
        self._max_gen = max_gen
        self._max_change = max_change

    def get_parameters(self):
        raise NotImplementedError

    def initialize(self):
        raise NotImplementedError

    def generate_chromosome_info(self, chromosome):
        raise NotImplementedError

    def evaluate(self, chromosome):
        raise NotImplementedError

    def select_and_reproduce(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError


class AssignmentsGA(GeneticAlgorithm):
    """
    Defines a Genetic Algorithm to assign nurses to patients
    """

    def __init__(self, rooms, distance_matrix, patients, pool_size=500, p_cross=0.8, p_muta=0.05, elitism=True,
                 endogamy=False, max_gen=5000, max_change=1500, nurses=4, w_loads=1, w_dist=1):
        super().__init__(pool_size, p_cross, p_muta, elitism, endogamy, max_gen, max_change)
        self._rooms = rooms
        self._distance_matrix = distance_matrix
        self._patients = patients
        self._nurses = nurses
        self._w_loads = w_loads
        self._w_dist = w_dist
        self._reproduce_p = [self._p_cross, 1 - self._p_cross]
        self._mutate_p = [self._p_muta, 1 - self._p_muta]
        self._verbose = True

    def get_parameters(self):
        """
        Retrieves the values of different parameters from the genetic algorithm
        :return: a string with all the required information
        """
        pool_size = str("Pool size: {0}\n".format(self._pool_size))
        p_cross = str("Crossover probability: {0}\n".format(self._p_cross))
        p_muta = str("Mutation probability: {0}\n".format(self._p_muta))
        elitism = str("Elitism enabled: {0}\n".format(self._elitism))
        endogamy = str("Endogamy enabled: {0}\n".format(self._endogamy))
        max_gen = str("Maximum generations: {0}\n".format(self._max_gen))
        max_change = str("Maximum generations without change: {0}\n".format(self._max_change))
        w_loads = str("Loads coeficient: {0}\n".format(self._w_loads))
        w_dist = str("Distances coeficient: {0}\n".format(self._w_dist))
        return pool_size + p_cross + p_muta + elitism + endogamy + max_gen + max_change + w_loads + w_dist

    def initialize(self):
        """
        Creates the first pool to evaluate
        :return: a random pool of chromosmes
        """
        base_population = [i for i in range(len(self._patients))]

        for nurse in range(-1, -self._nurses, -1):
            base_population.append(nurse)

        starting_pool = pool.Pool()

        for population in range(self._pool_size):
            gene_sequence = sample(base_population, len(base_population))
            starting_chromosome = chroms.AssignmentChromosome(gene_sequence)
            starting_chromosome.set_evaluation(self.evaluate(starting_chromosome))
            starting_pool.add_chromosome(starting_chromosome)

        starting_pool.calculate_enter_p()

        self._pool = starting_pool

    def generate_chromosome_info(self, chromosome):
        """
        Generates the information about a chromosome so it's readable
        :param chromosome: the crhomosome that needs its information polled
        :return: a string that contains the fitting score as well as the assignments
        """
        header = "Chromosome with score {0}\n".format(chromosome.get_evaluation())
        assigns = chromosome.get_assignments()
        loads = self.transform_assignments_loads(assigns)
        distances = self.calculate_distances(assigns)

        def generate_assignment_info(x, y, z):
            x = str(sorted(x)).strip("[]")
            y = str(sum(y))
            z = str(z)
            return "Assignment: {0}. Load score: {1}. Distance score: {2}".format(x, y, z)

        info_list = list(map(generate_assignment_info, assigns, loads, distances))

        info = ""

        if len(info_list) > 0:
            for element in info_list:
                info = info + str(element) + "\n"

        return header + info

    def evaluate(self, chromosome):
        """
        Evaluates a chromosome taking into account the fitting function
        :param chromosome: chromosome to be evaluated
        :return: the value of the fitting function for that chromosome
        """
        nurse_assignments = chromosome.get_assignments()
        loads = self.calculate_load(nurse_assignments)
        distance = sum(self.calculate_distances(nurse_assignments))

        evaluation = (self._w_loads * loads) + (self._w_dist * distance)
        return evaluation

    def calculate_load(self, assignments):
        """
        Calculates the load of care from the different assignments of patients
        :param assignments: sets of assignments
        :return: the standard variance of the cost of cares
        """

        loads = self.transform_assignments_loads(assignments)
        loads = [sum(load) for load in loads]
        mean = sum(loads) / len(loads)
        sq_differences = [(x - mean) ** 2 for x in loads]
        variance = sum(sq_differences) / len(sq_differences)

        return sqrt(variance)

    def transform_assignments_loads(self, assignments):
        return [[int(self._patients[index]) for index in assign] for assign in assignments]

    def calculate_distances(self, assignments):
        """
        Calculates the total distance that all nurses will need to connect all patients
        :param assignments: sets of patient assignments to nurses
        :return: the sum of the vaule of the minimal trees
        """
        nurse_distances = []
        for assign in assignments:
            nurse_distances.append(self.get_pst(assign))
        return nurse_distances

    def get_pst(self, assign):
        """
        Calculates the cost of the primordial tree
        :param assign: an assignment of patients for a nurse
        :return: the minimal cost of the tree that joins all patients
        """
        if len(assign) != 0:
            included_nodes = {}
            priority_queue = []
            heapq.heappush(priority_queue, (0, 0, 0))
            while len(included_nodes) != len(assign):
                minimum_node = heapq.heappop(priority_queue)
                distance = minimum_node[0]
                visiting_node = minimum_node[1]
                parent = minimum_node[2]
                if visiting_node not in included_nodes:
                    included_nodes[visiting_node] = (distance, parent)
                    for vertex in range(len(assign)):
                        if vertex not in included_nodes:
                            heapq.heappush(priority_queue,
                                           (self._distance_matrix[assign[vertex]][assign[visiting_node]],
                                            vertex, visiting_node))
            distance = 0
            for values in included_nodes.values():
                distance += values[0]
            return distance
        else:
            return 0

    def select_and_reproduce(self):
        """
        Creates a new pool of children chromosome from the actual pool
        :return: A new pool with child chromosomes
        """

        if self._elitism:
            candidates_list = choices(self._pool.get_pool_list(), weights=self._pool.get_enter_p(), k=self._pool_size - 1)
        else:
            candidates_list = choices(self._pool.get_pool_list(), weights=self._pool.get_enter_p(), k=self._pool_size)

        new_pool = pool.Pool()

        for candidate in candidates_list:
            parent = cp.copy(candidate)
            reproduce = choices(["cross", "no_change"], weights=self._reproduce_p)
            if reproduce[0] == "cross":
                if not self._endogamy:
                    partner = cp.copy(choices(self._pool.get_pool_list(), weights=self._pool.get_enter_p(), k=1)[0])
                else:
                    partner = self._pool.get_best_chromosome()
                crossover_type = choice(["ordered", "simple", "double"])
                if crossover_type == "ordered":
                    parent.ordered_crossover(partner)
                elif crossover_type == "simple":
                    parent.simple_crossover(partner)
                elif crossover_type == "double":
                    parent.double_crossover(partner)
            elif reproduce[0] == "no_change":
                pass

            mutate = choices(["muta", "no_change"], weights=self._mutate_p)
            if mutate[0] == "muta":
                parent.mutate()
            elif mutate[0] == "no_change":
                pass

            parent.set_evaluation(self.evaluate(parent))
            new_pool.add_chromosome(parent)

        if self._elitism:
            new_pool.add_chromosome(cp.copy(self._pool.get_best_chromosome()))

        new_pool.calculate_enter_p()

        self._pool = new_pool

    def run(self):
        start = time.time()
        base_file_name = create_file_name()

        parameters = open(base_file_name + "_param_asig", "w")
        parameters.write(self.get_parameters())
        parameters.close()

        chromosomes = open(base_file_name + "_chrom_asig", "w")

        log = open(base_file_name + "_log_asig", "w")
        log.write("MIN, MAX, SUM, MEAN, MEDIAN\n")

        self.initialize()

        generation_number = 0
        generations_no_changes = 0

        best_chromosome = cp.copy(self._pool.get_best_chromosome())
        chromosomes.write(self.generate_chromosome_info(best_chromosome))
        min_evaluation = best_chromosome.get_evaluation()

        log.write(str(self._pool.get_stats()).strip("()") + "\n")

        while (generation_number < self._max_gen) and (generations_no_changes < self._max_change):
            self.select_and_reproduce()
            if min_evaluation <= min(self._pool.get_evaluations()):
                generations_no_changes += 1
            else:
                min_evaluation = min(self._pool.get_evaluations())
                best_chromosome = cp.copy(self._pool.get_best_chromosome())
                chromosomes.write(self.generate_chromosome_info(best_chromosome))
                if self._verbose:
                    print("Better chromosome found at generation {0} with score {1}".format(
                        generation_number, best_chromosome.get_evaluation()))
                generations_no_changes = 0
            generation_number += 1

            log.write(str(self._pool.get_stats()).strip("()") + "\n")

            if generation_number % 100 == 0 and self._verbose:
                stats = self._pool.get_stats()
                print("Generation {0}: {1} min score, {2} max score, {3} sum score, {4} mean score, {5} median score"
                      .format(generation_number, stats[0], stats[1], stats[2], stats[3], stats[4]))
        end = time.time()
        if self._verbose:
                print("The best chromosome is {0} with score {1} after {2} generations in {3} seconds"
                      .format(best_chromosome.get_gene_sequence(), min_evaluation, generation_number, round(end-start)))

        chromosomes.close()
        log.close()
        return best_chromosome


class UbicationGA(GeneticAlgorithm):
    """
    Defines a Genetic Algorithm to explore better ubications for patients
    """
    def __init__(self, rooms, distance_matrix, patients, pool_size=100, p_cross=0.8, p_muta=0.05, elitism=True,
                 endogamy=True, max_gen=5000, max_change=1500, nurses=4, w_loads=1, w_dist=1):
        super().__init__(pool_size, p_cross, p_muta, elitism, endogamy, max_gen, max_change)
        self._rooms = rooms
        self._distance_matrix = distance_matrix
        self._patients = patients
        self._nurses = nurses
        self._w_loads = w_loads
        self._w_dist = w_dist
        self._verbose = True

    def get_parameters(self):
        pool_size = str("Pool size: {0}\n".format(self._pool_size))
        p_cross = str("Crossover probability: {0}\n".format(self._p_cross))
        p_muta = str("Mutation probability: {0}\n".format(self._p_muta))
        elitism = str("Elitism enabled: {0}\n".format(self._elitism))
        endogamy = str("Endogamy enabled: {0}\n".format(self._endogamy))
        max_gen = str("Maximum generations: {0}\n".format(self._max_gen))
        max_change = str("Maximum generations without change: {0}\n".format(self._max_change))
        w_loads = str("Loads coeficient: {0}\n".format(self._w_loads))
        w_dist = str("Distances coeficient: {0}\n".format(self._w_dist))
        return pool_size + p_cross + p_muta + elitism + endogamy + max_gen + max_change + w_loads + w_dist

    def initialize(self):
        """
        Creates the first pool of chromosomes
        :return: a pool with chromosomes with a swapping of patients
        """

        starting_pool = pool.Pool()

        processing_pool = mp.Pool(processes=mp.cpu_count())

        initial_chromosome_list = [chroms.UbicationChromosome([i for i in range(len(self._patients))])
                                   for individual in range(len(self._patients)-1)]

        starting_chromosome_list = processing_pool.map(self.mutate_cromosome, [c for c in initial_chromosome_list])

        processing_pool.close()
        processing_pool.join()

        for chromosome in starting_chromosome_list:
            starting_pool.add_chromosome(chromosome)

        starting_chromosome = chroms.UbicationChromosome([i for i in range(len(self._patients))])
        assignment_chromosome = self.evaluate(starting_chromosome)
        starting_chromosome.set_assignment_chromosome(assignment_chromosome)
        starting_pool.add_chromosome(starting_chromosome)

        starting_pool.calculate_enter_p()

        self._pool = starting_pool

    def mutate_cromosome(self, chromosome):
        chromosome.mutate()
        assignment_chromosome = self.evaluate(chromosome)
        chromosome.set_assignment_chromosome(assignment_chromosome)
        return chromosome

    def generate_chromosome_info(self, chromosome):
        header = "Chromosome with score {0}\n".format(chromosome.get_evaluation())
        info = "{0} \n \n".format(chromosome.get_gene_sequence())
        return header + info

    def evaluate(self, chromosome):
        """
        Calculates the best Assignment Chromosome for an ubication Chromosome
        :param chromosome: the patient distribution from where to calculate its value
        :return: the value of the chromosome
        """
        new_patient_list = self.map_loads(chromosome.get_gene_sequence())
        assignment_chromosome = AssignmentsGA(self._rooms, self._distance_matrix, new_patient_list,
                                              nurses=self._nurses).run()
        return assignment_chromosome

    def map_loads(self, gene_sequence):
        loads = [self._patients[index] for index in gene_sequence]
        return loads

    def select_and_reproduce(self):
        """
        Creates a new pool of child chromosomes
        :return: a new pool from where to keep calculating
        """
        if self._elitism:
            candidates_list = choices(self._pool.get_pool_list(), weights=self._pool.get_enter_p(), k=self._pool_size - 1)
        else:
            candidates_list = choices(self._pool.get_pool_list(), weights=self._pool.get_enter_p(), k=self._pool_size)

        new_pool = pool.Pool()

        new_list = []

        for candidate in candidates_list:
            parent = cp.copy(candidate)
            new_list.append(parent)

        processing_pool = mp.Pool(processes=mp.cpu_count())

        processing_pool.map(self.mutate_cromosome, [c for c in new_list])

        processing_pool.close()
        processing_pool.join()

        for child in new_list:
            new_pool.add_chromosome(child)

        if self._elitism:
            new_pool.add_chromosome(cp.copy(self._pool.get_best_chromosome()))

        new_pool.calculate_enter_p()

        self._pool = new_pool

    def run(self):
        start = time.time()
        base_file_name = create_file_name()

        parameters = open(base_file_name + "_param_ubi", "w")
        parameters.write(self.get_parameters())
        parameters.close()

        chromosomes = open(base_file_name + "_chrom_ubi", "w")

        log = open(base_file_name + "_log_ubi", "w")
        log.write("MIN, MAX, SUM, MEAN, MEDIAN\n")

        self.initialize()

        generation_number = 0
        generations_no_changes = 0

        best_chromosome = cp.copy(self._pool.get_best_chromosome())
        chromosomes.write(self.generate_chromosome_info(best_chromosome))
        min_evaluation = best_chromosome.get_evaluation()

        log.write(str(self._pool.get_stats()).strip("()") + "\n")

        while (generation_number < self._max_gen) and (generations_no_changes < self._max_change):
            self.select_and_reproduce()
            if min_evaluation <= min(self._pool.get_evaluations()):
                generations_no_changes += 1
            else:
                min_evaluation = min(self._pool.get_evaluations())
                best_chromosome = cp.copy(self._pool.get_best_chromosome())
                chromosomes.write(self.generate_chromosome_info(best_chromosome))
                if self._verbose:
                    print("Better chromosome found at generation {0} with score {1}".format(
                        generation_number, best_chromosome.get_evaluation()))
                generations_no_changes = 0
            generation_number += 1

            log.write(str(self._pool.get_stats()).strip("()") + "\n")

            if generation_number % 100 == 0 and self._verbose:
                stats = self._pool.get_stats()
                print("Generation {0}: {1} min score, {2} max score, {3} sum score, {4} mean score, {5} median score"
                      .format(generation_number, stats[0], stats[1], stats[2], stats[3], stats[4]))
        end = time.time()
        if self._verbose:
            print("The best chromosome is {0} with score {1} after {2} generations in {3} seconds"
                  .format(best_chromosome.get_gene_sequence(), min_evaluation, generation_number, round(end-start)))

        chromosomes.close()
        log.close()
        return best_chromosome
