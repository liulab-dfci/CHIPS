#!/usr/bin/env python

from genson import SchemaBuilder
import json

###############################################################################
# Constants
###############################################################################
_int = 777
_float = 7.0
_string = "TGBTG"
_int_arr = [_int,_int,_int]
_float_arr = [_float, _float, _float]
_file_path = "/the/only/way"

###############################################################################
# SAMPLE level definition
###############################################################################
alignment = {'total_reads': _int,
             'mapped_reads': _int,
             'dedup_reads': _int,
             'gc_content': _float_arr, #NOTE: could be floats
             'quality_score': _float_arr } #NOTE: could be floats


contamination = {
  'contamination_percent' : _float_arr #NOTE: could be floats
}

pbc = {'Nd_count': _int,
       'N1_count': _int,
       'pbc': _float } # Note: could be floats


sample = {'id': _string,
          'raw_file1': _file_path, #bucket path, either fastq or bam
          'raw_file2': _file_path, #bucket path, in case of PE
          'dedup_bam_file': _file_path,
          'alignment': alignment,
          'contamination': contamination,
          'pbc': pbc
          }

###############################################################################
# Immunoprep level definition
###############################################################################
frip = {'reads_in_peaks': _int,
        'total_reads': _int,
        'frip': _float}

peaks_qc = {'total_peaks' : _int,
         'fc_10': _int,
         'fc_20': _int,
         'dhs': _float,
         'promoter': _int,
         'exon': _int,
         'intron': _int,
         'intergenic': _int }

targets = {'peak_targets': _file_path }


conservations = {'x_axis': _int_arr,
                 'y_axis': _float_arr }


ip_results = {'peaks': _file_path,
              'peaks_qc': peaks_qc,
              'frip': frip,
              'targets': targets,
              'conservations': conservations}

###############################################################################
# RUN level definition
###############################################################################


run = {'id': _string,
       'ip': sample,
       'control': sample,
       'ip_results': ip_results}


###############################################################################
# END RUN level definition
###############################################################################

builder = SchemaBuilder(schema_uri="http://json-schema.org/draft-07/schema#")
builder.add_object(run)
print(builder.to_json(indent=3))

js = builder.to_json(indent=3)
with open("chips_schema.json", "w") as outfile:
  outfile.write(js)





