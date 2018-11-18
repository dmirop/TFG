#!/usr/bin/python

from random import sample, choice, randrange


def mutation(chromosome):
    mutation_points = sample(range(len(chromosome)), 2)

    mutated_gene_a = chromosome[mutation_points[0]]
    mutated_gene_b = chromosome[mutation_points[1]]

    mutated_chromosome = chromosome.copy()

    mutated_chromosome[mutation_points[0]] = mutated_gene_b
    mutated_chromosome[mutation_points[1]] = mutated_gene_a

    return mutated_chromosome


def simple_crossover(parent_a, parent_b, c_point=None, c_type=None):
    if c_point is None:
        c_point = randrange(len(parent_a))

    if c_type is None:
        c_type = choice(["start", "end"])

    start_parent_a = parent_a[:c_point]
    end_parent_a = parent_a[c_point:]

    start_parent_b = parent_b[:c_point]
    end_parent_b = parent_b[c_point:]

    if c_type == "start":
        missing_genes = [gene for gene in end_parent_a if gene not in end_parent_b]
        conflict_index = [index for index in range(len(end_parent_b)) if end_parent_b[index] in start_parent_a]

        def assign_gene(x, y):
            end_parent_b[x] = y

        list(map(assign_gene, conflict_index, missing_genes))

        crossover_chromosome = [*start_parent_a, *end_parent_b]

    elif c_type == "end":
        missing_genes = [gene for gene in start_parent_a if gene not in start_parent_b]
        conflict_index = [index for index in range(len(start_parent_b)) if start_parent_b[index] in end_parent_a]

        def assign_gene(x, y):
            start_parent_b[x] = y

        list(map(assign_gene, conflict_index, missing_genes))

        crossover_chromosome = [*start_parent_b, *end_parent_a]

    else:
        raise AttributeError("Not a valid crossover type: {0}".format(c_type))

    return crossover_chromosome


def double_crossover(parent_a, parent_b, crossover_points=None, crossover_type=None):
    if crossover_points is None:
        crossover_points = []

    if len(crossover_points) == 0:
        crossover_points = sorted(sample(range(len(parent_a)), 2))
    elif len(crossover_points) != 2:
        raise AttributeError("Not a valid number of crossover points: crossover points length {0}"
                             .format(len(crossover_points)))

    if crossover_type is None:
        crossover_type = choice(["middle", "borders"])

    start_parent_a = parent_a[:crossover_points[0]]
    middle_parent_a = parent_a[crossover_points[0]:crossover_points[1]]
    end_parent_a = parent_a[crossover_points[1]:]

    start_parent_b = parent_b[:crossover_points[0]]
    middle_parent_b = parent_b[crossover_points[0]:crossover_points[1]]
    end_parent_b = parent_b[crossover_points[1]:]

    borders_parent_a = [*start_parent_a, *end_parent_a]
    borders_parent_b = [*start_parent_b, *end_parent_b]

    if crossover_type == "borders":

        missing_genes = [gene for gene in middle_parent_a if gene not in middle_parent_b]
        conflict_index = [index for index in range(len(middle_parent_b)) if
                          middle_parent_b[index] in borders_parent_a]

        def assign_gene(x, y):
            middle_parent_b[x] = y

        list(map(assign_gene, conflict_index, missing_genes))

        crossover_chromosome = [*start_parent_a, *middle_parent_b, *end_parent_a]

    else:

        missing_genes = [gene for gene in borders_parent_a if gene not in borders_parent_b]
        conflict_index = [index for index in range(len(borders_parent_b)) if borders_parent_b[index] in middle_parent_a]

        def assign_gene(x, y):
            borders_parent_b[x] = y

        list(map(assign_gene, conflict_index, missing_genes))

        crossover_chromosome = [*borders_parent_b[:crossover_points[0]], *middle_parent_a,
                                *borders_parent_b[crossover_points[1]-len(parent_a):]]

    return crossover_chromosome


def ordered_crossover(parent_a, parent_b):
    nurse_index_a = [i for i, x in enumerate(parent_a) if x < 0]

    part_keep = choice(range(len(nurse_index_a)))

    if part_keep == 0:
        return simple_crossover(parent_a, parent_b, nurse_index_a[0], c_type="start")
    elif part_keep == len(nurse_index_a):
        return simple_crossover(parent_a, parent_b, nurse_index_a[-1], c_type="end")
    else:
        return double_crossover(parent_a, parent_b, [nurse_index_a[part_keep - 1],
                                                     nurse_index_a[part_keep] + 1], "middle")
