"""Tools"""
from .. import dicty
from Orange.data import ContinuousVariable, StringVariable, TimeVariable, Domain, Table


CallBack = dicty.CallBack
transformValues = dicty.transformValues
averageAttributes = dicty.averageAttributes
example_tables = dicty.example_tables


def transpose_table(table):
    """
    Transpose the rows and columns of the table.

    Args:
        table: Data in :obj:`Orange.data.Table`

    Returns:
         Transposed :obj:`Orange.data.Table`. (Genes as columns)
    """
    attrs = table.domain.attributes
    attr = [ContinuousVariable.make(ex['Gene'].value) for ex in table]
    #  Set metas
    new_metas = [StringVariable.make(name) if name is not 'Time' else TimeVariable.make(name)
                 for name in sorted(table.domain.variables[0].attributes.keys())]
    domain = Domain(attr, metas=new_metas)
    meta_values = [[exp.attributes[var.name] for var in domain.metas] for exp in attrs]

    return Table(domain, table.X.transpose(), metas=meta_values)


def to_orange_table(expres_dict, anno, view_model, transform=False, transpose=False,
                        joinreplicates=False, callback=lambda: None):


    """ Transform downloaded data to :obj:`Orange.data.Table`.

    Args:
        expres_dict (dict): Dictionary of downloaded expressions. Key represents a gene and
                            values are time-points.

        anno (list): Annotations for each sample.
                     E.g. [['Dd_Dp_genome_biology', 'K.a.', 'wildtype', 'Prespore',
                         'NC4', 16, 2, '2015-04-20T13:28:07.836000']]

        view_model (list): Annotations to include as attributes in the table.
                           E.g. ['Experiment', 'Growth', 'Genotype', 'Treatment', 'Strain', 'Time', 'Replicate']

        transform (bool): Logarithmic (base 2) transformation. Default: False
        transpose (bool): Represent genes al columns.
        joinreplicates (bool): Join replicates by median.

    Returns:
        Data in :obj:`Orange.data.Table`.

    """

    variables = []
    table = []
    cbc = CallBack(1, callback, callbacks=10)
    for i in anno:
        var = ContinuousVariable(i[0])  # name of the variable is the same as sample name
        for n in list(zip(view_model, i))[1:]:
            var.attributes[n[0]] = str(n[1])
        variables.append(var)
    cbc()

    meta_attr = StringVariable.make('Gene')
    domain = Domain(variables, metas=[meta_attr])

    for gen, time in expres_dict.items():
        time.append(gen)
        table.append(time)

    orange_table = Table(domain, table)
    cbc.end()

    cbc = CallBack(3, callback, callbacks=10)
    if transform:
        transformValues(orange_table, fn=transform)  # in place transform
        cbc()

    if joinreplicates:
        orange_table = dicty.join_replicates(orange_table,
                                             namefn="name", avg=dicty.median,
                                             fnshow=lambda x: " | ".join(map(str, x)),
                                             ignorenames=["Replicate"])
        cbc()

    # sort output columns
    view_model_dict = {name: name for name in view_model}
    if 'Time' in view_model_dict:
        def sorting_key(attr):
            sortOrder = ['Time']
            atts = attr.attributes
            return tuple([int(atts.get(view_model_dict[name])) for name in sortOrder])


        attributes = sorted(orange_table.domain.attributes, key=sorting_key)

        domain = Domain(attributes, metas=[meta_attr])
        orange_table = Table.from_table(domain, orange_table)
        cbc()

    if transpose:
        orange_table = transpose_table(orange_table)

    cbc.end()

    return orange_table
