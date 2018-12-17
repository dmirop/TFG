#!/usr/bin/python


class Pool:
    def __init__(self):
        self._pool_list = []
        self._best_chromosome = None
        self._enter_p = None

    def add_chromosome(self, chromosome):
        """
        Method to add a new chromosome to the pool
        :param chromosome: chromosome to add
        """
        self._pool_list.append(chromosome)
        if self._best_chromosome is None:
            self._best_chromosome = chromosome
        elif chromosome.get_evaluation() < self._best_chromosome.get_evaluation():
            self._best_chromosome = chromosome

    def get_best_chromosome(self):
        return self._best_chromosome

    def get_evaluations(self):
        return [chromosome.get_evaluation() for chromosome in self._pool_list]

    def min_evaluation(self):
        return self._best_chromosome.get_evaluation()

    def max_evaluation(self):
        return max([chromosome.get_evaluation() for chromosome in self._pool_list])

    def get_stats(self):
        """
        Method to retrieve the statistics of the pool
        :return: a tuple of statistics
        """
        evaluations = self.get_evaluations()
        return min(evaluations), max(evaluations), sum(evaluations), sum(evaluations)/len(evaluations), \
            sorted(evaluations)[round(len(evaluations)/2)]

    def get_pool_list(self):
        return self._pool_list

    def set_pool(self, pool_list):
        self._pool_list = pool_list

    def calculate_enter_p(self):
        """
        Calculates the probability to be chosen for the next generation pool
        :return: a list with the probabilities
        """
        evaluations = self.get_evaluations()
        sum_evaluations = sum(evaluations)

        if sum_evaluations > 0:
            converting_factor = max(evaluations) + min(evaluations)
            self._enter_p = [round((converting_factor - evaluation) / sum_evaluations, 4) for evaluation in evaluations]
        else:
            self._enter_p = evaluations

    def get_enter_p(self):
        return self._enter_p
