from abc import ABC, abstractmethod


# TODO Complete this file
class FileWriter(ABC):

    @classmethod
    def __init__(cls, filename):
        cls.filename = filename

    @abstractmethod
    def _write(self, filename):
        pass


class PrintGmt(FileWriter):
    """
    This class imlpements the writing of a GMT file on a file
    """
    def __init__(self, filename):
        super().__init__(filename)

    def _write(self, filename):
        pass
