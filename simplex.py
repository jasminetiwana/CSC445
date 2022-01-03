# Jasmine Tiwana V00850594
#!/bin/env python3
from typing import List, Dict, Tuple
import sys
import numpy as np

def main():
    array1 = []
    for line in sys.stdin:
        line.strip()
        list = line.split()
        array1.append(list)
    table = generate_table(array1)
    variables = []
    for i in range(1, ((np.shape(table)[1]-1) + np.shape(table)[0]-1)):
        variables += [i]
    status = simplexmethod(variables, table)
    if (status == -2):
        status = generate_aux(table)
    if (status == -1):
        print('Unbounded')
    if (status == -2):
        print('Infeasible')

def generate_table(array1):
    constraints = array1[1:]
    list = []
    for constraint in constraints:
        list.append(constraint.pop())
    for constraint in constraints:
        for i in range (0, len(constraint)):
            constraint[i] = float(constraint[i])
            constraint[i] *= (-1)
    obj = array1[0]
    objective = np.array([None, 0])
    objective = np.hstack(([None], [0], obj))
    basis = np.array([0] * len(constraints))
    for i in range(0, len(basis)):
        basis[i] = len(obj) + i + 1
    basis2 = np.hstack((np.transpose([basis]), np.transpose([list]), constraints))
    table = np.vstack((objective, basis2))
    table = np.array(table, dtype=float)

    return table

def enteringvar(function):
    curr = 2
    largeind = -1
    largest = 0
    for num in function[0][2:]:
        if num > 0:
            if num > largest:
                largest = num
                largeind = curr
        curr += 1
    return largeind

def leavingvar(table, entering):
    leaving = -1
    min = 99999
    for i in range(1, len(table)):
        c = table[i][1]
        if c < 0:
            return -2
        if (table[i][entering] <= -0.00000000000000001) and (table[i][1] > 0):
            m = table[i][entering]
            c = table[i][1]
            if c < 0:
                return -2 #infeasible
            val = abs(c/m)
            if val < min:
                min = val
                leaving = i
    return leaving

def generate_aux(table):
    obj = table[0]
    omega = np.ones((np.shape(table)[0],1))
    aux = np.append(table, omega, axis=1)
    aux[0] *= 0
    aux[0][-1] = -1
    entering = len(aux[0])-1
    min = 1
    for i in range(1, np.shape(aux)[0]):
        if (aux[i][1] < min):
            min = aux[i][1]
            leaving = i
    aux_variables = []
    for i in range(1, ((np.shape(aux)[1]-1) + np.shape(aux)[0]-1)):
        aux_variables.append(i)
    for i in range(1, np.shape(table)[0]):
        aux[i][0] += 1
    pivot = aux[leaving][entering]
    aux[leaving][0] = aux_variables[entering - 2]
    aux_variables[entering-2], aux_variables[leaving + (len(aux[0])-2) - 1] = aux_variables[leaving + (len(aux[0])-2) - 1], aux_variables[entering-2]
    for i in range(1, len(aux[0])):
        if (i == entering):
            aux[leaving][i] = 1
        else:
            aux[leaving][i] *= -1
    for i in range(0, len(aux)):
        if (i != leaving):
            if aux[i][entering] != 0:
                mult = aux[leaving] * aux[i][entering]
                for j in range(1, len(aux[0])):
                    if (j == entering):
                        aux[i][entering] /= pivot
                    else:
                        aux[i][j] += mult[j]
    optimal = False
    while (True):
        #pivoting
        entering = enteringvar(aux)
        leaving = leavingvar(aux, entering)
        if (leaving == -2):
            return -2 #infeasible
        if (entering == -1):
            optimal = True
        if optimal == True:
            break
        if leaving == -1:
            return -1
        pivot = aux[leaving][entering]
        aux[leaving][0] = aux_variables[entering - 2]
        aux_variables[entering-2], aux_variables[leaving + (len(aux[0])-2) - 1] = aux_variables[leaving + (len(aux[0])-2) - 1], aux_variables[entering-2]
        for i in range(1, len(aux[0])):
            if (i == entering):
                aux[leaving][i] = 1/pivot
            else:
                aux[leaving][i] = aux[leaving][i]/ abs(pivot)
        for i in range(0, len(aux)):
            if (i != leaving):
                if aux[i][entering] != 0:
                    mult = aux[leaving] * aux[i][entering]
                    for j in range(1, len(aux[0])):
                        if (j == entering):
                            aux[i][entering] /= pivot
                        else:
                            aux[i][j] += mult[j]
    if (optimal == True):
        if (aux[0][1] != 0):
            return (-2) #infeasible
    for item in aux_variables:
        if item == (len(aux[0])-2):
            index = aux_variables.index(item)
    aux_variables.remove(len(table[0])-1)
    table = np.delete(aux, index+2, 1)
    for i in range(1, np.shape(table)[0]):
        if table[i][0] > len(table[0]):
            table[i][0] -= 1
    obj_variables = len(obj)-2
    obj_index = []
    for i in range(1, obj_variables+1):
        for item in aux_variables:
            if i == item:
                index = aux_variables.index(item)
                obj_index.append(index)
    obj_function = []
    for item in obj_index:
        if item >= (len(table[0])-2):
            for i in range(1, np.shape(table)[0]):
                if table[i][0] == aux_variables[item]:
                    function = table[i] * obj[obj_index.index(item) + 2]
                    obj_function.append(function)
    obj_function = np.array(obj_function)
    obj_function2 = np.sum(obj_function, 0)
    table[0] = obj_function2
    table[0][0] = None
    simplexmethod(aux_variables, table)

def simplexmethod(variables, table):
    optimal = False
    while (True):
        #pivoting
        entering = enteringvar(table)
        leaving = leavingvar(table, entering)
        if (leaving == -2):
            break #infeasible
        if (leaving == -1):
            break
        if (entering == -1):
            optimal = True
        if optimal == True:
            break
        pivot = table[leaving][entering]
        table[leaving][0] = variables[entering - 2]
        variables[entering-2], variables[leaving + (len(table[0])-2) - 1] = variables[leaving + (len(table[0])-2) - 1], variables[entering-2]
        for i in range(1, len(table[0])):
            if (i == entering):
                table[leaving][i] = 1/pivot
            else:
                table[leaving][i] = table[leaving][i]/ abs(pivot)
        for i in range(0, len(table)):
            if (i != leaving):
                if table[i][entering] != 0:
                    mult = table[leaving] * table[i][entering]
                    for j in range(1, len(table[0])):
                        if (j == entering):
                            table[i][entering] /= pivot
                        else:
                            table[i][j] += mult[j]
    if (optimal == True):
        print('optimal')
        print(table[0][1])
        for i in range(1, len(table[0])-1):
            for item in variables:
                if i == item:
                    if variables.index(item) >= (len(table[0])-2):
                        print(table[variables.index(item)-(len(table[0])-3)][1], end = ' ')
                    else:
                        print(0, end = ' ')
        print()
        return
    if (entering == -1 and leaving == -1):
        return -2
    if (leaving == -2):
        return -2
    if (leaving == -1):
        return -1

main()
