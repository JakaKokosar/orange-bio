import Orange
import numpy


unknown = NAN = float("nan")
isunknown = numpy.isnan


def create_domain(at, cl, metas):
    return Orange.data.Domain(at, cl, metas=metas)


def create_table(domain, X, Y, metas):
    classvar = domain.class_var
    metaatts = domain.metas

    if Y:
        Y = numpy.array([[classvar.to_val(row)] for row in Y],
                        dtype=float)
    if metas:
        metas = numpy.array([[c.to_val(v) for c, v in zip(metaatts, row)]
                             for row in metas],
                            dtype=object)
    data = Orange.data.Table(domain, numpy.asarray(X), Y=Y, metas=metas)

    return data

