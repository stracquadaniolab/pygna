import pygna.reading_class as rc


# TODO Check and refactor this
def __load_geneset(filename, setname):
    """Loads a geneset from file
    """
    return rc.ReadGmt(filename).get_geneset(setname)
