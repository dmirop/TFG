#!/usr/bin/python

import sys
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
            exit(0)
        return True
    except AssertionError:
        print("error: the number of arguments is not correct \n")
        print("Try 'TFG --help' for more information")
        exit(1)


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
    return abs(point_a[0] - point_b[0]) + abs(point_a[1] - point_b[1])


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


def main():
    global rooms
    global distances
    global patients

    if check_arguments():

        retrieve_rooms()

        calculate_room_distances()

        retrieve_patients()

    genetic_algorithm = ga.AssignmentsGA(rooms, distances, patients, pool_size=50, nurses=2)
    best_chromosome = genetic_algorithm.run()


if __name__ == "__main__":
    main()
