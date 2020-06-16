import pandas as pd

runTable = pd.read_csv("ChIPRunTable.csv", sep = ",")

runTable['aggregate_column'] = runTable[['BioSample', 'library_type']].agg('_'.join, axis = 1)
runTable.to_csv("ChIPRunTable.csv", sep = ",", index = False)