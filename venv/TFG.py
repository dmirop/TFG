#!/usr/bin/python

import sys
import random
from math import sqrt
import heapq
from itertools import accumulate


def check_arguments():
    """Asserts the right number of arguments

    :return: Indicates if the assertion is right
    """
    try:
        assert 1 < len(sys.argv) < 4
        if len(sys.argv) == 2:
            assert sys.argv[1] == "--help"
            print("Usage: TFG ROOMS PATIENTS")
            print("Executes a genetic algorithm taking into account the ROOMS and PATIENTS")
            exit(0)
        return True
    except AssertionError:
        print("error: the number of arguments is not correct \n")
        print("Try 'TFG --help' for more information")


def retrieve_rooms():
    """Tries to retrieve the rooms from file

    :return: A list with all the coordinates identifying each room
    """
    try:
        f = open(sys.argv[1], 'r')
    except OSError:
        print('Cannot open', sys.argv[1])
        print("Try 'TFG --help' for more information")
    else:
        rooms = []
        for t in f.read().split():
            a, b = t.strip('()').split(',')
            rooms.append((int(a), int(b)))
        f.close()
        return rooms


def calculate_room_distances(rooms):
    """Computes the distance between rooms

    :param rooms: a list of room coordinates
    :return: a list of distance lists
    """
    distances_list = []
    for from_room in rooms:
        distance_room_list = []
        for to_room in rooms:
            distance_room_list.append(manhattan_distance(from_room, to_room))
        distances_list.append(distance_room_list)

    return distances_list


def manhattan_distance(point_a, point_b):
    """Gives the Manhattan distance of two coordinates

    :param point_a: coordinates of point A
    :param point_b: coordinates of point B
    :return: the Manhattan distance
    """
    return abs(point_a[0] - point_a[0]) + abs(point_a[1] - point_b[1])


def retrieve_patients():
    """Tries to retrieve the patients from a file

    :return: a list of patients
    """
    try:
        f = open(sys.argv[2], 'r')
    except OSError:
        print('Cannot open', sys.argv[2])
        print("Try 'TFG --help' for more information")
    else:
        patients = f.read().splitlines()
        f.close()
        return patients


def initialize(population_size, num_patients, nurses):
    """
    Creates de first generation of chromosomes
    :param population_size: size of the population
    :param num_patients: number of patients
    :param nurses: number of nurses in each turn
    :return: a pool for the first generation
    """
    base_population = [i for i in range(num_patients)]

    for nurse in range(-1, -nurses, -1):
        base_population.append(nurse)

    pool = []
    for population in range(population_size):
        pool.append(random.sample(base_population, len(base_population)))
    return pool


def evaluate(pool_list, distance_matrix, patient_list, coeficient_stdv, coeficient_distance):
    evaluations = []
    for chromosome in pool_list:
        nurse_assignments = retrieve_assignments(chromosome)
        stdv = calculate_stdv(nurse_assignments, patient_list)
        distance = calculate_distances(nurse_assignments, distance_matrix)
        evaluations.append((coeficient_stdv * stdv) + (coeficient_distance * distance))
    return evaluations


def retrieve_assignments(chromosome):
    assigns = []
    nurse_assign = []
    for gene in chromosome:
        if gene >= 0:
            nurse_assign.append(gene)
        else:
            assigns.append(nurse_assign)
            nurse_assign = []
    assigns.append(nurse_assign)

    return assigns


def calculate_stdv(assignments, patient_list):
    n = len(assignments)
    loads = []
    for assign in assignments:
        load = 0
        if len(assign) != 0:
            for patient in assign:
                load += int(patient_list[patient])
        loads.append(load)
    mean = sum(loads) / n
    sq_differences = [(x - mean) ** 2 for x in loads]
    variance = sum(sq_differences) / n
    return sqrt(variance)


def calculate_distances(assignments, distances):
    nurse_distances = []
    for assign in assignments:
        nurse_distances.append(get_PST(assign, distances))
    return sum(nurse_distances)


def get_PST(assign, distances):
    """
    Calculate Prim's Spanning Tree Algorithm
    :param assign: the assignment from where to calculate the Spanning Tree
    :param distances: distance matrix where to look up the necessary distances
    :return: the path value of the assignment
    """
    INF = 99999
    INF = 0

    if len(assign) != 0:
        included_nodes = {}
        priority_queue = []
        heapq.heappush(priority_queue, (0, 0, 0))
        while len(included_nodes) != len(assign):
            minimum_node = heapq.heappop(priority_queue)
            distance = minimum_node[0]
            visiting_node = minimum_node[1]
            parent = minimum_node[2]
            if visiting_node not in [*included_nodes]:
                included_nodes[visiting_node] = (distance, parent)
                for vertex in range(len(assign)):
                    if vertex not in [*included_nodes]:
                        heapq.heappush(priority_queue,
                                       (distances[assign[vertex]][assign[visiting_node]], vertex, visiting_node))
        distance = 0
        for values in included_nodes.values():
            distance += values[0]
        return distance
    else:
        return INF


def select_and_reproduce(pool, mutation_rate, crossover_rate, evaluations):
    enter_p = [(max(evaluations)+min(evaluations))-eval for eval in evaluations]
    enter_p = [round(eval/sum(enter_p),4) for eval in enter_p]
    enter_p = list(accumulate(enter_p))
    enter_p[-1] = 1

    best_chromosome = pool[evaluations.index(min(evaluations))]

    change_p = [mutation_rate, mutation_rate+crossover_rate,1]

    new_pool = []

    for i in range(len(pool)-1):
        rndm = random.random()
        for j in range(len(enter_p)):
            if rndm <= enter_p[j]:
                new_pool.insert(i, pool[j])
                break

    for i in range(len(new_pool)):
        rndm = random.random()
        for j in range(len(change_p)):
            if rndm <= change_p[j]:
                if j == 0:
                    new_pool[i] = mutation(new_pool[i])
                elif j == 1:
                    crossover_type = random.choice(["oc", "sc", "dc"])
                    if crossover_type == "oc":
                        new_pool[i] = ordered_crossover(new_pool[i], best_chromosome)
                    elif crossover_type == "sc":
                        new_pool[i] = simple_crossover(new_pool[i], best_chromosome)
                    elif crossover_type == "dc":
                        new_pool[i] = double_crossover(new_pool[i], best_chromosome)
                break

    new_pool.append(best_chromosome)

    return new_pool

def mutation(chromosome):
    mutation_points = random.sample(range(len(chromosome)),2)

    mutated_gene_A = chromosome[mutation_points[0]]
    mutated_gene_B = chromosome[mutation_points[1]]

    mutated_chromosome = chromosome.copy()

    mutated_chromosome[mutation_points[0]] = mutated_gene_B
    mutated_chromosome[mutation_points[1]] = mutated_gene_A

    return mutated_chromosome


def ordered_crossover(parent_a, parent_b):
    nurse_index_a = [i for i, x in enumerate(parent_a) if x < 0]

    part_keep = random.choice(range(len(nurse_index_a)))

    if part_keep == 0:
        return simple_crossover(parent_a, parent_b, nurse_index_a[0])
    elif part_keep == len(nurse_index_a):
        return simple_crossover(parent_a, parent_b, nurse_index_a[-1])
    else:
        return double_crossover(parent_a, parent_b, [nurse_index_a[part_keep - 1], nurse_index_a[part_keep] + 1], "borders")


def simple_crossover(parent_a, parent_b, crossover_point = None):
    if crossover_point == None:
        crossover_point = random.randrange(len(parent_a))

    end_parent_b = parent_b[crossover_point:]
    missing_genes = [gene for gene in parent_a[crossover_point:] if gene not in end_parent_b]
    conflict_index = [index for index in range(len(end_parent_b))if end_parent_b[index] in parent_a[:crossover_point]]

    def assign_gene(x, y):
        end_parent_b[x] = y

    list(map(assign_gene, conflict_index, missing_genes))

    crossover_chromosome = [*parent_a[:crossover_point],*end_parent_b]

    return crossover_chromosome

def double_crossover(parent_a, parent_b, crossover_points=[], crossover_type= None):
    if len(crossover_points) == 0:
        crossover_points = sorted(random.sample(range(len(parent_a)), 2))

    if crossover_type == None:
        crossover_type = random.choice(["middle", "borders"])

    if crossover_type == "middle":

        middle_parent_b = parent_b[crossover_points[0]:crossover_points[1]]

        missing_genes = [gene for gene in parent_a[crossover_points[0]:crossover_points[1]] if
                         gene not in middle_parent_b]
        conflict_index = [index for index in range(len(middle_parent_b)) if
                          middle_parent_b[index] in parent_a[:crossover_points[0]] or middle_parent_b[
                              index] in parent_a[crossover_points[1]:]]

        def assign_gene(x, y):
            middle_parent_b[x] = y

        list(map(assign_gene, conflict_index, missing_genes))

        crossover_chromosome = [*parent_a[:crossover_points[0]], *middle_parent_b, *parent_a[crossover_points[1]:]]

    else:
        start_parent_a = parent_a[:crossover_points[0]]
        middle_parent_a = parent_a[crossover_points[0]:crossover_points[1]]
        end_parent_a = parent_a[crossover_points[1]:]


        start_parent_b = parent_b[:crossover_points[0]]
        middle_parent_b = parent_b[crossover_points[0]:crossover_points[1]]
        end_parent_b = parent_b[crossover_points[1]:]

        border_genes_a = [*start_parent_a, *end_parent_a]
        border_genes_b = [*start_parent_b, *end_parent_b]

        missing_genes = [gene for gene in border_genes_a if gene not in border_genes_b]
        conflict_index = [index for index in range(len(border_genes_b)) if border_genes_b[index] in middle_parent_a]

        def assign_gene(x, y):
            border_genes_b[x] = y

        list(map(assign_gene, conflict_index, missing_genes))

        crossover_chromosome = [*border_genes_b[:crossover_points[0]], *middle_parent_a, *border_genes_b[crossover_points[1]-len(parent_a):]]

    return crossover_chromosome


def main():
    if check_arguments():
        rooms = retrieve_rooms()
        distances = calculate_room_distances(rooms)
        patients = retrieve_patients()
    num_nurses = 4
    pop_size = 100
    crossover_prob = 0.05
    mutation_prob = 0.9
    max_generations = 50000
    gen_no_change = max_generations/100

    weight_stdv = 1
    weight_distance = 1

    pool = initialize(pop_size, len(patients), num_nurses)
    evaluations = evaluate(pool, distances, patients, weight_stdv, weight_distance)
    generation_number = 0
    generations_no_changes = 0
    min_evaluation = min(evaluations)
    best_chromosome = pool[evaluations.index(min_evaluation)].copy()

    while (generation_number < max_generations) and (generations_no_changes < gen_no_change):
        pool = select_and_reproduce(pool, mutation_prob, crossover_prob, evaluations)
        evaluations = evaluate(pool, distances, patients, weight_stdv, weight_distance)
        if min_evaluation <= min(evaluations):
            generations_no_changes += 1
        else:
            min_evaluation = min(evaluations)
            best_chromosome = pool[evaluations.index(min_evaluation)].copy()
            nurse_assignments = retrieve_assignments(best_chromosome)
            distance = calculate_distances(nurse_assignments, distances)
            stdv = calculate_stdv(nurse_assignments, patients)
            print(
                "Better chromosome found at generation {0} with stdv {1} and distance {2}".format(
                    generation_number, stdv, distance))
            generations_no_changes = 0
        generation_number += 1

        if generation_number%100 == 0:
            print("Generation {0}: {1} min score, {2} max score, {3} sum score, {4} mean score, {5} median score".format(generation_number, min(evaluations), max(evaluations), sum(evaluations), sum(evaluations)/len(evaluations), sorted(evaluations)[round(len(evaluations)/2)]))

    print("The best chromosome is {0} with score {1}".format(best_chromosome, min_evaluation))

if __name__ == "__main__":
    main()