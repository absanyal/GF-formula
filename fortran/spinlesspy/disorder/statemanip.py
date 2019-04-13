# Thu Feb 7 11:33:37 IST 2019
eta = 0.1
t = complex(1.0, 0.0)
u_int = 8.0  # Change in BOTH files


class state:

    def __init__(self, s=0, ns=0, alphas=0, betas=0, totalu=0):
        self.s = s
        self.ns = ns
        self.alphas = alphas
        self.betas = betas
        self.totalu = totalu


def kdelta(x, y):
    if (x == y):
        return 1
    else:
        return 0


def factorial(n):
    p = 1
    while (n > 1):
        p = p * n
        n -= 1
    return p


def setstate(astate, vs, vns, valphas, vbetas, vtotalu):
    astate.s = vs
    astate.ns = vns
    astate.alphas = valphas
    astate.betas = vbetas
    astate.totalu = vtotalu


def getstate(astate):
    return [astate.s, astate.ns, astate.alphas, astate.betas]


def checkvalidity(s1):
    cv = 1
    if (s1.s < 1 or s1.ns < 0):
        cv = 0
    if (s1.s == 0 and s1.ns == 0 and
            s1.alphas == 0 and s1.betas == 0):
        cv = 0
    if (s1.s - 1 < s1.ns - s1.alphas):
        cv = 0
    if (s1.alphas > s1.ns):
        cv = 0
    return cv


def getstatesize(s1):
    gs = 1
    vs = s1.s
    vns = s1.ns
    va = s1.alphas
    vb = s1.betas
    if (checkvalidity(s1) == 1):
        if (va == 0 and vb == 0):
            gs = \
                factorial(vs - 1) / (factorial(vns) * factorial(vs - 1 - vns))
        if (va == 0 and vb == 1):
            gs = \
                factorial(vs - 1) / (factorial(vns) * factorial(vs - 1 - vns))
        if (va == 1 and vb == 0):
            gs = \
                factorial(vs - 1) / (factorial(vns - 1) * factorial(vs - vns))
        if (va == 1 and vb == 1):
            gs = \
                factorial(vs - 1) / (factorial(vns - 1) * factorial(vs - vns))
    else:
        gs = 0

    return int(gs)


def relegate0(s1):
    r0 = state()
    sm1 = s1.s - 1
    betasm1 = s1.alphas
    alphasm1 = 0
    nsm1 = s1.ns - betasm1
    tusm1 = s1.totalu
    setstate(r0, sm1, nsm1, alphasm1, betasm1, tusm1)
    return r0


def relegate1(s1):
    r1 = state()
    sm1 = s1.s - 1
    betasm1 = s1.alphas
    alphasm1 = 1
    nsm1 = s1.ns - betasm1
    tusm1 = s1.totalu
    if (betasm1 == 1):
        tusm1 += u_int
    setstate(r1, sm1, nsm1, alphasm1, betasm1, tusm1)
    return r1


def getlstate(s1):
    if (s1.alphas == 0 and s1.betas == 0):
        k = "a"
    if (s1.alphas == 0 and s1.betas == 1):
        k = "b"
    if (s1.alphas == 1 and s1.betas == 0):
        k = "c"
    if (s1.alphas == 1 and s1.betas == 1):
        k = "d"
    return str(int(s1.s)) + k + str(int(s1.ns))


def splitstates(fname, s1):
    if (checkvalidity(s1) == 1):
        # print(getlstate(s1), '\t', getstatesize(s1))
        p = str(getlstate(s1)) + '\t' \
            + str(getstatesize(s1)) + '\t'\
            + str(s1.totalu) + '\t' + \
            '\n'
        fname.write(p)
        s1s0 = relegate0(s1)
        splitstates(fname, s1s0)
        s1s1 = relegate1(s1)
        splitstates(fname, s1s1)


def splitstatesatlevel(fname, s1, level):
    if (checkvalidity(s1) == 1):
        if (s1.s == level):
            p = str(getlstate(s1)) + '\t' \
                + str(getstatesize(s1)) + '\t'\
                + str(s1.totalu) + '\n'
            fname.write(p)
        if (s1.s > level):
            s1s0 = relegate0(s1)
            splitstatesatlevel(fname, s1s0, level)
            s1s1 = relegate1(s1)
            splitstatesatlevel(fname, s1s1, level)


def splitstatesatlevelraw(fname, s1, level):
    if (checkvalidity(s1) == 1):
        if (s1.s == level):
            fname.write(
                str(getstate(s1)[0]) + '\t' + str(getstate(s1)[1]) + '\t' +
                str(getstate(s1)[2]) + '\t' + str(getstate(s1)[3]) + '\t'
                + str(getstatesize(s1))+'\n'
            )
        if (s1.s > level):
            s1s0 = relegate0(s1)
            splitstatesatlevelraw(fname, s1s0, level)
            s1s1 = relegate1(s1)
            splitstatesatlevelraw(fname, s1s1, level)


def getletter(s1):
    if (s1.alphas == 0 and s1.betas == 0):
        return 1
    if (s1.alphas == 0 and s1.betas == 1):
        return 2
    if (s1.alphas == 1 and s1.betas == 0):
        return 3
    if (s1.alphas == 1 and s1.betas == 1):
        return 4

########################################################


def z(w):
    import numpy as np
    return (w + np.complex(0, eta))


def G(H, w):
    import numpy as np
    return np.linalg.inv(z(w)*np.eye(len(H), dtype=np.complex) - H)


def tmm(A, B, C):
    import numpy as np
    return np.dot(A, np.dot(B, C))

########################################################


def connecting_u(s1, s2):
    import numpy as np
    if (checkvalidity(s1) == 1 and checkvalidity(s2) == 1):
        u = \
            np.zeros((getstatesize(s1), getstatesize(s2)),
                     dtype=np.complex)
        s1s0 = relegate0(s1)
        startindex = getstatesize(s1s0)
        i = startindex
        j = 0
        while(i < getstatesize(s1)):
            u[i, j] = np.complex(1, 0)
            i += 1
            j += 1
    return u
