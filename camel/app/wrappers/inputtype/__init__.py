from camel.app.wrappers.inputtype.fastahelper import FastaHelper
from camel.app.wrappers.inputtype.illuminahelper import IlluminaHelper
from camel.app.wrappers.inputtype.onthelper import ONTHelper

helper_by_input_type = {
    'fasta': FastaHelper,
    'illumina': IlluminaHelper,
    'ont': ONTHelper
}
