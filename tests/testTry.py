import pandas as pd

try:
    pd.load('templates/pulsar_parameters')
except IOError:
    print "DNE"