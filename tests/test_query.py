""" These tests are obsolete (we know the prototype works) """

def test_queries(self):
    """Test all queries and see if they correspond to the polynomial evaluation.
    """

    m, d, q = self.m, self.d, self.q

    error = False
    fast_poly = fast_callable(self.polynomial, vars=self.variables)
    for i in range(q**(m-1)):
        # resolving data
        point = self.interpolating_set[i]
        
        server_no, idx= self.resolve_query(point)
        symbol = self.carry_out_query(server_no, idx, 0)
        
        if fast_poly(*point) != symbol:
            error = True
    if not error:
        print('Success')


def test_queries(self):
    """Test all queries and see if they correspond to the polynomial evaluation.
    """

    comm = Communicator(self.mode, self.m, self.d, self.q, self.s)

    m, d, q, s = self.m, self.d, self.q, self.s

    poly_indices_dict = comp_indices_dict(self.encoder.polynomial, m)

    partials = []
    for i in constrained_bounded_partitions(s, m):
        partial = vector(i)
        partial.set_immutable()
        partials.append(partial)

    error = False
    # generate all derivatives of the polynomial
    derivs = []
    for chunk, partial in enumerate(partials):
            deriv = hasse_derivative(self.encoder.polynomial.monomials(), poly_indices_dict, self.encoder.polynomial.coefficients(), partial, self.encoder.variables, m)
            fast_deriv = fast_callable(deriv, vars=self.encoder.variables)
            derivs.append(fast_deriv)
    # evaluate all derivatives of the polynomial
    for i in range(q**(m-1)):
        point = self.encoder.interpolating_set[i]
        server_no, idices = comm.resolve_query(point)
        symbols = comm.carry_out_query(server_no, idices, 0)
        evals = []
        for chunk, partial in enumerate(partials):
            # resolving data
            evals.append(derivs[chunk](*point))
                
        if evals != symbols:
            error = True
            print('ERROR')

    if not error:
        print('Successful querying')
    else:
        print('Unsuccessful querying')
