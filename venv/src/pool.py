#!/usr/bin/python


class Pool:
    def __init__(self):
        self._pool_list = []
        self._best_chromosome = None

    def add_chromosome(self, chromosome):
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
        evaluations = self.get_evaluations()
        return min(evaluations), max(evaluations), sum(evaluations), sum(evaluations)/len(evaluations), sorted(evaluations)[round(len(evaluations)/2)]

    def get_pool_list(self):
        return self._pool_list

    def set_pool(self, pool_list):
        self._pool_list = pool_list
