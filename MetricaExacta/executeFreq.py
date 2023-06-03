import subprocess

genomas = [10, 12, 19, 29, 33, 35, 38, 40, 41, 45, 46, 47, 52, 56, 62,69,71,72,73,76,79,82,83,89,92,96,98]
for g in genomas:
    for k in [10, 15]:
        bashCommand = "./freq_reales GCF_00{g} {k}".format(g=g, k=k)
        try:
            subprocess.run(bashCommand, shell=True, check=True, capture_output=True)
            string = bashCommand.split(' ')
            print("Genomas: {string1}, Ejecucion k={k}.-".format(string1=string[1], k=k))
        except subprocess.CalledProcessError as e:
                string=bashCommand.split(' ')
                print("Genomas: {string1}, Ejecucion k={k}.-".format(string1=string[1], k=k))
