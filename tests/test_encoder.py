""" These tests are obsolete (we know the prototype works) """

def test_interpolation(self):

    m, d, q, s, k = self.m, self.d, self.q, self.s, self.k
    sigma, k_out = self.sigma, self.k_out
    Q = self.polynomial

    message_seq = self.translate_message()
    message_segments = [list() for _ in range(sigma)]
    for i in range(k_out):
        message_segments[i % sigma].append(message_seq[i])
    
    for i in range(sigma):
        message_segments[i] = vector(message_segments[i])

    partials = []
    for i in constrained_bounded_partitions(s, m):
        partial = vector(i)
        partial.set_immutable()
        partials.append(partial)

    Q_indices_dict = comp_indices_dict(Q, m)

    error = False
    for idx, i in enumerate(partials):
        deriv = hasse_derivative(Q.monomials(), Q_indices_dict, Q.coefficients(), i, self.variables, m)
        fast_deriv = fast_callable(deriv, vars=self.variables)
        evals = []
        for p in self.interpolating_set:
            eval_pt = fast_deriv(*p)
            evals.append(eval_pt)
        evals = vector(evals)
        if evals != message_segments[idx]:
            error = True
    if not error:
        print('Successful interpolation')
    else:
        print('Unsuccessful interpolation')