################
import sys, os
from fisher_py.data import Device, TrayShape
from fisher_py.raw_file_reader import RawFileReaderAdapter
################

def load_raw(filename):
    rawFile = RawFileReaderAdapter.file_factory(
                                                os.getcwd()+"/src/../../raw_files/"+filename
                                                )
    rawFile.select_instrument(Device.MS, 1)
    return rawFile

