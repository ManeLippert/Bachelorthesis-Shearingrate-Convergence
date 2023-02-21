import pandas as pd
import ast

df = pd.read_csv('data.csv', delimiter=';')

start = ast.literal_eval(df['evo_start'][16])

print(start)
print(df)