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
        return rooms
        f.close()


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


def manhattan_distance(A, B):
    """Gives the Manhattan distance of two coordinates

    :param A: coordinates of point A
    :param B: coordinates of point B
    :return: the Manhattan distance
    """
    return abs(A[0] - B[0]) + abs(A[1] - B[1])


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
        return f.read().splitlines()
        f.close()


def initialize(popSize, num_patients):
    """
    Creates de first generation of chromosomes
    :param popSize: size of the population
    :param num_patients: number of patients
    :return: a pool for the first generation
    """
    base_population = [i for i in range(num_patients)]

    for nurse in range(2):
        base_population.append("N")

    pool = []
    for population in range(popSize):
        pool.append(random.sample(base_population, len(base_population)))
    return pool


def evaluate(pool_list, distance_matrix, patient_list):
    evaluations = []
    for chromosome in pool_list:
        nurse_assignments = retrieve_assignments(chromosome)
        stdv = calculate_stdv(nurse_assignments, patient_list)
        #distance = calculate_distances(nurse_assignments, distance_matrix)
        distance = 0
        evaluations.append(stdv + 2*distance)
    return evaluations


def retrieve_assignments(chromosome):
    assigns = []
    nurse_assign = []
    for gene in chromosome:
        if gene != "N":
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

    if len(assign) != 0:
        included_nodes = {}
        priority_queue = []
        heapq.heappush(priority_queue, (0, 0, 0))
        while len(included_nodes) != len(assign):
            minimum_node = heapq.heappop(priority_queue)
            distance = minimum_node[0]
            visiting_node = minimum_node[1]
            parent = minimum_node[2]
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

    enter_p = [sum(evaluations)-eval for eval in evaluations]
    enter_p = [round(eval/sum(enter_p),4) for eval in enter_p]
    enter_p = list(accumulate(enter_p))
    enter_p[-1] = 1

    best_chromosome = pool[evaluations.index(min(evaluations))]

    change_p = [mutation_rate, mutation_rate+crossover_rate,1]

    new_pool = []

    for i in range(len(pool)):
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
                    new_pool[i] = crossover(new_pool[i], best_chromosome)
                break

    return new_pool

def crossover(chromosome_a, chromosome_b):
    nurse_index_a = get_nurse_index(chromosome_a)

    kept_assign = random.randrange(len(nurse_index_a) + 1)

    if kept_assign == 0:
        protected_range = [0, nurse_index_a[kept_assign] - 1]
    elif kept_assign == len(nurse_index_a):
        protected_range = [nurse_index_a[kept_assign - 1] + 1, len(chromosome_a) - 1]
    else:
        protected_range = [nurse_index_a[kept_assign - 1] + 1, nurse_index_a[kept_assign] - 1]

    crossover_chromosome = []
    copy_chromosome_b = chromosome_b.copy()

    for i in range(len(chromosome_a)):
        if chromosome_a[i] == "N":
            crossover_chromosome.insert(i, "N")
        elif protected_range[0] <= i <= protected_range[1]:
            crossover_chromosome.insert(i, chromosome_a[i])
        else:
            while len(crossover_chromosome) != i + 1:
                candidate_gene = copy_chromosome_b.pop(0)
                if candidate_gene != "N":
                    if candidate_gene not in crossover_chromosome:
                        if candidate_gene not in chromosome_a[protected_range[0]:protected_range[1] + 1]:
                            crossover_chromosome.insert(i, candidate_gene)
    return crossover_chromosome


def get_nurse_index(chromosome):
    return [i for i, x in enumerate(chromosome) if x == "N"]


def mutation(chromosome):
    gene_index = []
    while len(gene_index) != 2:
        index = random.randrange(len(chromosome))
        if index not in gene_index:
            gene_index.append(index)

    mutated_gene_A = chromosome[gene_index[0]]
    mutated_gene_B = chromosome[gene_index[1]]

    chromosome[gene_index[0]] = mutated_gene_B
    chromosome[gene_index[1]] = mutated_gene_A

    return chromosome


def main():
    if check_arguments():
        rooms = retrieve_rooms()
        distances = calculate_room_distances(rooms)
        patients = retrieve_patients()

    pop_size = 100
    crossover_prob = 0.2
    mutation_prob = 0.2
    max_generations = 10000
    gen_no_change = 250

    pool = initialize(pop_size, len(patients))
    evaluations = evaluate(pool, distances, patients);
    generation_number = 0
    generations_no_changes = 0
    min_evaluation = min(evaluations)
    best_chromosome = pool[evaluations.index(min_evaluation)].copy()

    while (generation_number < max_generations) and (generations_no_changes < gen_no_change):
        pool = select_and_reproduce(pool, mutation_prob, crossover_prob, evaluations)
        evaluations = evaluate(pool, distances, patients);
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

        if generation_number%10 == 0:
            print("Generation {0}: {1} min score".format(generation_number, min(evaluations)))

    print("The best chromosome is {0} with score {1}".format(best_chromosome, min_evaluation))


if __name__ == "__main__":
    main()
