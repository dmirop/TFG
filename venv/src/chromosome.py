#!/usr/bin/python

from random import sample, randrange, choice


class Chromosome:
    def __init__(self, gene_sequence):
        self._gene_sequence = gene_sequence
        self._evaluation = 0

    def get_evaluation(self):
        return self._evaluation

    def get_gene_sequence(self):
        return self._gene_sequence

    def set_evaluation(self, evaluation):
        self._evaluation = evaluation

    def set_gene_sequence(self, gene_sequence):
        self._gene_sequence = gene_sequence

    def mutate(self):
        """
        Method to mutate a chromosome
        :return: BaseException if there has been an error
        """
        previous_gene_seq = [*self._gene_sequence]
        mutation_points = sample(range(len(self._gene_sequence)), 2)

        mutated_gene_a = self._gene_sequence[mutation_points[0]]
        mutated_gene_b = self._gene_sequence[mutation_points[1]]

        self._gene_sequence[mutation_points[0]] = mutated_gene_b
        self._gene_sequence[mutation_points[1]] = mutated_gene_a

        if sorted(self._gene_sequence) != sorted(previous_gene_seq):
            raise BaseException("The mutation has produced a non valid chromosome")

    def simple_crossover(self, other, c_point=None, c_type=None):
        """
        Method that produces a new child with simple crossover
        :param other: the partner used for the crossover
        :param c_point: crossover point
        :param c_type: type of crossover
        :return: a new child with new gene sequence
        """
        if c_point is None:
            c_point = randrange(len(self._gene_sequence))

        if c_type is None:
            c_type = choice(["start", "end"])

        previous_gene_seq = [*self._gene_sequence]

        start_self = self._gene_sequence[:c_point]
        end_self = self._gene_sequence[c_point:]

        start_other = other.get_gene_sequence()[:c_point]
        end_other = other.get_gene_sequence()[c_point:]

        if c_type == "start":
            missing_genes = [gene for gene in end_self if gene not in end_other]
            conflict_index = [index for index in range(len(end_other)) if end_other[index] in start_self]

            def assign_gene(x, y):
                end_other[x] = y

            list(map(assign_gene, conflict_index, missing_genes))

            self.set_gene_sequence([*start_self, *end_other])
            if sorted(self._gene_sequence) != sorted(previous_gene_seq):
                raise BaseException("The simple start crossover has produced a non valid chromosome")

        elif c_type == "end":
            missing_genes = [gene for gene in start_self if gene not in start_other]
            conflict_index = [index for index in range(len(start_other)) if start_other[index] in end_self]

            def assign_gene(x, y):
                start_other[x] = y

            list(map(assign_gene, conflict_index, missing_genes))

            self.set_gene_sequence([*start_other, *end_self])
            if sorted(self._gene_sequence) != sorted(previous_gene_seq):
                raise BaseException("The simple end crossover has produced a non valid chromosome")

        else:
            raise AttributeError("Not a valid crossover type: {0}".format(c_type))

    def double_crossover(self, other, crossover_points=None, crossover_type=None):
        """
        Method that produces a new child with double crossover
        :param other: the partner used for the crossover
        :param c_point: crossover points
        :param c_type: type of crossover
        :return: a new child with new gene sequence
        """
        if crossover_points is None:
            crossover_points = []

        if len(crossover_points) == 0:
            crossover_points = sorted(sample(range(len(self._gene_sequence)), 2))
        elif len(crossover_points) != 2:
            raise AttributeError("Not a valid number of crossover points: crossover points length {0}"
                                 .format(len(crossover_points)))

        if crossover_type is None:
            crossover_type = choice(["middle", "borders"])

        previous_gene_seq = [*self._gene_sequence]

        start_parent_a = self._gene_sequence[:crossover_points[0]]
        middle_parent_a = self._gene_sequence[crossover_points[0]:crossover_points[1]]
        end_parent_a = self._gene_sequence[crossover_points[1]:]

        start_parent_b = other.get_gene_sequence()[:crossover_points[0]]
        middle_parent_b = other.get_gene_sequence()[crossover_points[0]:crossover_points[1]]
        end_parent_b = other.get_gene_sequence()[crossover_points[1]:]

        borders_parent_a = [*start_parent_a, *end_parent_a]
        borders_parent_b = [*start_parent_b, *end_parent_b]

        if crossover_type == "borders":

            missing_genes = [gene for gene in middle_parent_a if gene not in middle_parent_b]
            conflict_index = [index for index in range(len(middle_parent_b)) if
                              middle_parent_b[index] in borders_parent_a]

            def assign_gene(x, y):
                middle_parent_b[x] = y

            list(map(assign_gene, conflict_index, missing_genes))

            self.set_gene_sequence([*start_parent_a, *middle_parent_b, *end_parent_a])
            if sorted(self._gene_sequence) != sorted(previous_gene_seq):
                raise BaseException("The double borders crossover has produced a non valid chromosome")

        else:
            missing_genes = [gene for gene in borders_parent_a if gene not in borders_parent_b]
            conflict_index = [index for index in range(len(borders_parent_b)) if
                              borders_parent_b[index] in middle_parent_a]

            def assign_gene(x, y):
                borders_parent_b[x] = y

            list(map(assign_gene, conflict_index, missing_genes))

            self.set_gene_sequence([*borders_parent_b[:crossover_points[0]], *middle_parent_a,
                                    *borders_parent_b[crossover_points[0]:]])
            if sorted(self._gene_sequence) != sorted(previous_gene_seq):
                raise BaseException("The double middle crossover has produced a non valid chromosome")

    def print(self):
        return "Chromosome with score {0}\n".format(self.get_evaluation())


class AssignmentChromosome(Chromosome):
    """
    Class to define an assignment Chromosome
    """
    def ordered_crossover(self, other):
        """
        Creates a new child with ordered crossover
        :param other: the other parent that participates
        :return: a new child with new gene sequence
        """
        nurse_index = self.get_nurse_indices()

        part_keep = choice(range(len(nurse_index)))

        if part_keep == 0:
            self.simple_crossover(other, nurse_index[0], c_type="start")
        elif part_keep == len(nurse_index):
            self.simple_crossover(other, nurse_index[-1], c_type="end")
        else:
            self.double_crossover(other, [nurse_index[part_keep - 1],
                                          nurse_index[part_keep] + 1], "middle")

    def get_assignments(self):
        """
        Method to retrieve the nurse assignments
        :return: a list of assignments
        """
        assigns = []
        nurse_assign = []
        for gene in self._gene_sequence:
            if gene >= 0:
                nurse_assign.append(gene)
            else:
                assigns.append(nurse_assign)
                nurse_assign = []
        assigns.append(nurse_assign)

        return assigns

    def get_nurse_indices(self):
        return [i for i, x in enumerate(self._gene_sequence) if x < 0]

