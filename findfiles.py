import os

pattern = ".ipynb"

for dirpath, dirnames, filenames in os.walk("."):
    for filename in [f for f in filenames if f.endswith(pattern)]:
        print(dirpath + "/" + filename)