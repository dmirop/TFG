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


def ordered_crossover(parent_a, parent_b):
    nurse_index_a = [i for i, x in enumerate(parent_a) if x < 0]

    part_keep = choice(range(len(nurse_index_a)))

    if part_keep == 0:
        return simple_crossover(parent_a, parent_b, nurse_index_a[0])
    elif part_keep == len(nurse_index_a):
        return simple_crossover(parent_a, parent_b, nurse_index_a[-1])
    else:
        return double_crossover(parent_a, parent_b, [nurse_index_a[part_keep - 1],
                                                     nurse_index_a[part_keep] + 1], "borders")


def simple_crossover(parent_a, parent_b, crossover_point=None):
    if crossover_point is None:
        crossover_point = randrange(len(parent_a))

    end_parent_b = parent_b[crossover_point:]
    missing_genes = [gene for gene in parent_a[crossover_point:] if gene not in end_parent_b]
    conflict_index = [index for index in range(len(end_parent_b))if end_parent_b[index] in parent_a[:crossover_point]]

    def assign_gene(x, y):
        end_parent_b[x] = y

    list(map(assign_gene, conflict_index, missing_genes))

    crossover_chromosome = [*parent_a[:crossover_point], *end_parent_b]

    return crossover_chromosome


def double_crossover(parent_a, parent_b, crossover_points=None, crossover_type=None):
    if crossover_points is None:
        crossover_points = []
    if len(crossover_points) == 0:
        crossover_points = sorted(sample(range(len(parent_a)), 2))

    if crossover_type is None:
        crossover_type = choice(["middle", "borders"])

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

        crossover_chromosome = [*border_genes_b[:crossover_points[0]], *middle_parent_a,
                                *border_genes_b[crossover_points[1]-len(parent_a):]]

    return crossover_chromosome
