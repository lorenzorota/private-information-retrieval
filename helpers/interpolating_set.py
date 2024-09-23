from sage.all import *
from helpers.math import *
from helpers.misc import *
import random
import concurrent.futures


class InterpolatingSet:
    """ Container for generating appropriate interpolating sets and Vandermonde matrices
    given some coding parameters, and for writing and reading both objects into and from files.

    """

    def __init__(self, code, disjoint='False'):

        if type(code).__name__ == 'ReedMullerEncoder' or type(code).__name__ == 'ReedMullerLocalDecoder':
            self.mode = 'ReedMuller'
        elif type(code).__name__ == 'MultiplicityEncoder' or type(code).__name__ == 'MultiplicityLocalDecoder':
            self.mode = 'Multiplicity'
            self.s = code.s
            self.k = code.k
        else:
            raise TypeError('Object is not an encoder')
        if type(code).__name__ == 'ReedMullerEncoder' or type(code).__name__ == 'MultiplicityEncoder':
            self.monomials = code.monomials
            self.variables = code.variables
        self.m = code.m
        self.d = code.d
        self.q = code.q
        self.k = code.k
        self.vandermondes = []
        self.interpolating_sets = []


    def comp_vandermonde_concurrent(self, obj):
        vm_matrix, points = obj
        vm_matrix = copy(vm_matrix)
        for pos in range(20):
            vm_matrix[pos] = vm_matrix[pos].subs(dict(zip(self.variables, points[pos])))
        return vm_matrix


    def comp_vandermonde(self, points=None, concurrency=False):
        """ Evaluates a Vandermonde matrix from a given monomial basis and list of points.
        
        Args:
            points: A set of points as m-tuples.
            concurrency (bool): Indicate whether the computation should be performed over multiple processes (not that effective).
        
        Returns:
            vm_matrix: The evaluated generalized Vandermonde matrix

        """
        
        if points is None:
            points = self.interpolating_set[-1]

        q = self.q
        k = self.k

        if concurrency:
            processes = []
            new_matr = []
            process_count = 4
            upp = 40
            low = 20

            partition_size = 1
            for i in range(low,upp+1):
                common_div = gcd(k, i)
                if common_div == 1:
                    continue
                else:
                    partition_size = common_div
                    break

            vm_matrix = matrix([self.monomials for _ in range(partition_size)])
            executor = concurrent.futures.ProcessPoolExecutor(process_count)
            for points in grouper(points, partition_size):
                a = executor.submit(self.comp_vandermonde_concurrent, (vm_matrix, points))
                processes.append(a)
            concurrent.futures.wait(processes)
            for proc in processes:
                new_matr += list(proc.result())
            vm_matrix = matrix(GF(q), new_matr)
        else:
            vm_matrix = matrix(GF(q), [self.monomials.subs(dict(zip(self.variables, point))) for point in points])
        return vm_matrix


    def find_set(self, total_number=1, max_tries=100):
        """ Find a random interpolating set
        
        Each set is of size k = (m+d) choose d
        Every point in the set is an m-tuple
        
        Warning:
        This cannot feasibly be used in real-time if q^m or k are very large
        Instead, one should already know of appropriate interpolating sets
        
        Args:
            max_tries (default: False): The maximum number of iterations allowed.

        """
        
        Vandermonde_0 = None
        subset_0 = None
        self.interpolating_sets = []
        self.vandermondes = []
        m, q, k = self.m, self.q, self.k


        for i in range(total_number):
            for j in range(max_tries):
                # heuristic (otherwise this is impossible)
                subset = []
                for _ in range(k):
                    repeat = False
                    while not repeat:
                        point = vector(GF(q), (to_GF(q, random.randrange(1,q)) for _ in range(m)))
                        point.set_immutable()
                        if point in subset or any(point in other_sets for other_sets in self.interpolating_sets):
                            repeat = False
                        else:
                            repeat = True
                            subset.append(point)
                subset = tuple(subset)
                Vandermonde = self.comp_vandermonde(subset)

                # this is a lot faster than checking the determinant apparently
                if Vandermonde.inverse().inverse() == Vandermonde:
                    Vandermonde_0 = Vandermonde
                    subset_0 = subset
                    print("Found (disjoint) interpolating set: S_{} after {} tries ".format(i, j+1))
                    break
                
            if subset_0 is None:
                print("Nothing found")
            else:
                self.interpolating_sets.append(subset_0)
                self.vandermondes.append(Vandermonde_0)


    def GF_to_int_tuple(self, lst, m=0):
        """ Converts the field elements of the interpolating set from GF(q) to an integer modulo q """
        
        q, k = self.q, self.k
        new_lst = []
        
        if m == 0:
            for elem in lst:
                new_lst.append(from_GF(q, elem))
        else:
            for i in range(k):
                tup = []
                for j in range(m):
                    tup.append(from_GF(q, lst[i][j]))
                new_lst.append(tuple(tup))
        
        return tuple(new_lst)


    def int_to_GF_tuple(self, lst, m=0):
        """ Converts the integers of the interpolating set to elements from GF(q) """
            
        q, k = self.q, self.k
        new_lst = []
        
        if m == 0:
            for elem in lst:
                new_lst.append(to_GF(q, elem))
        else:
            for i in range(k):
                tup = []
                for j in range(m):
                    tup.append(to_GF(q, lst[i][j]))
                new_lst.append(tuple(tup))

        return tuple(new_lst)


    def build_keys(self):
        """Saves the interpolating set as a key
        
        Returns:
            key (str): A hex representation of the interpolating set

        """

        keys = []
        for interpolating_set in self.interpolating_sets:

            key = self.GF_to_int_tuple(interpolating_set, self.m)
            # create hex dump
            payload = ''
            for tup in key:
                for coord in tup:
                    payload += '{:0>2X}'.format(coord)
            key = payload
            keys.append(payload)
        
        return keys


    def build_interpolating_sets_from_keys(self, keys):
        """Load the key and generate the corresponding interpolating set."""
        
        m, k = self.m, self.k
        self.interpolating_set = []

        for key in keys:
            # create interpolating set
            lst_of_tuples = []
            for i in range(k):
                tup = []
                for j in range(m*i, m*(i+1)):
                    tup.append(int(key[2*j:2*j+2],16))
                lst_of_tuples.append(tuple(tup))
            
            self.interpolating_set.append(self.int_to_GF_tuple(lst_of_tuples, self.m))

        return self.interpolating_set