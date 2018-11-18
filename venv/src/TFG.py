#!/usr/bin/python

import sys
import random
from math import sqrt
import heapq
from itertools import accumulate
import chromosome
import pool
import ga

rooms = []
distances = []
patients = []


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
            exit(1)
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
        exit(1)
    else:
        for t in f.read().split():
            a, b = t.strip('()').split(',')
            rooms.append((int(a), int(b)))
        f.close()


def calculate_room_distances():
    """Computes the distance between rooms

    :param rooms: a list of room coordinates
    :return: a list of distance lists
    """
    for from_room in rooms:
        distance_room_list = []
        for to_room in rooms:
            distance_room_list.append(manhattan_distance(from_room, to_room))
        distances.append(distance_room_list)


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
        exit(1)
    else:
        for patient in f.read().splitlines():
            patients.append(patient)
        f.close()


def initialize(population_size, nurses):
    """
    Creates the first generation of chromosomes
    :param population_size: size of the population
    :param num_patients: number of patients
    :param nurses: number of nurses in each turn
    :return: a pool for the first generation
    """
    base_population = [i for i in range(len(patients))]

    for nurse in range(-1, -nurses, -1):
        base_population.append(nurse)

    pool = []
    for population in range(population_size):
        pool.append(random.sample(base_population, len(base_population)))
    return pool


def evaluate(pool_list, coeficient_stdv, coeficient_distance):
    evaluations = []
    for chromosome in pool_list:
        nurse_assignments = retrieve_assignments(chromosome)
        stdv = calculate_stdv(nurse_assignments)
        distance = calculate_distances(nurse_assignments)
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


def calculate_stdv(assignments):

    loads = transform_assignments_loads(assignments)

    loads = [sum(load) for load in loads]

    mean = sum(loads)/len(loads)

    sq_differences = [(x - mean) ** 2 for x in loads]

    variance = sum(sq_differences)/len(sq_differences)

    return sqrt(variance)


def transform_assignments_loads(assignments):

    return [[int(patients[index]) for index in assign] for assign in assignments]


def calculate_distances(assignments):
    nurse_distances = []
    for assign in assignments:
        nurse_distances.append(get_pst(assign))
    return sum(nurse_distances)


def get_pst(assign):
    """
    Calculate Prim's Spanning Tree Algorithm
    :param assign: the assignment from where to calculate the Spanning Tree
    :param distances: distance matrix where to look up the necessary distances
    :return: the path value of the assignment
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
                                       (distances[assign[vertex]][assign[visiting_node]], vertex, visiting_node))
        distance = 0
        for values in included_nodes.values():
            distance += values[0]
        return distance
    else:
        return 0


def select_and_reproduce(pool, mutation_rate, crossover_rate, evaluations):
    best_chromosome = pool[evaluations.index(min(evaluations))]

    converting_factor = max(evaluations) + min(evaluations)
    sum_evaluations = sum(evaluations)
    enter_p = [round((converting_factor-evaluation)/sum_evaluations, 4) for evaluation in evaluations]
    enter_p = list(accumulate(enter_p))

    change_p = [mutation_rate, crossover_rate+mutation_rate, 1]

    new_pool = random.choices(pool, cum_weights=enter_p, k=len(pool)-1)

    for index in range(len(new_pool)):
        change = random.choices(["mut", "cross", "no_change"], cum_weights=change_p)

        if change[0] == "mut":
            new_pool[index] = chromosome_utils.mutation(new_pool[index])

        elif change[0] == "cross":
            crossover_type = random.choice(["ordered", "simple", "double"])
            if crossover_type == "ordered":
                new_pool[index] = chromosome_utils.ordered_crossover(new_pool[index], best_chromosome)
            elif crossover_type == "simple":
                new_pool[index] = chromosome_utils.simple_crossover(new_pool[index], best_chromosome)
            elif crossover_type == "double":
                new_pool[index] = chromosome_utils.double_crossover(new_pool[index], best_chromosome)

        elif change[0] == "no_change":
            pass

    new_pool.append(best_chromosome)

    return new_pool


def main():
    global rooms
    global distances
    global patients

    if check_arguments():

        retrieve_rooms()

        calculate_room_distances()

        retrieve_patients()

    random.seed(0)

    genetic_algorithm = ga.AssignmentsGA(rooms, distances, patients)
    genetic_algorithm.run()

if __name__ == "__main__":
    main()
