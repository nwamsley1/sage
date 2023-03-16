################
import sys, os
from fisher_py.data import Device, TrayShape
from fisher_py.raw_file_reader import RawFileReaderAdapter
################

def load_raw(filename):
    rawFile = RawFileReaderAdapter.file_factory(
                                                #os.getcwd()+"/src/../../raw_files/"+filename
                                                #"/Users/n.t.wamsley/Projects/SAGE_TESTING/raw_files/" + filename
                                                #filename
                                                #"/Users/n.t.wamsley/Projects/SAGE_TESTING/MA4358_FFPE_HPVpos_01_071522.raw"
                                                filename
                                                )
    rawFile.select_instrument(Device.MS, 1)
    return rawFile

