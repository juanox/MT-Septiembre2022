import subprocess

genomas = [14, 15]
for g in genomas:
    for k in [10,15]:
        for p in [11, 12, 13]:
            bashCommand = "./entropy GCF_00{g} {k} {p} 3 2 11 2".format(g=g, k=k, p=p)
            try:
                subprocess.run(bashCommand, shell=True, check=True, capture_output=True)
                string = bashCommand.split(' ')
                print("Genomas: {string1}, Ejecucion k={k} p={p} .-".format(string1=string[1],  g=g,
                                                                                                                        k=k,
                                                                                                                        p=p
                                                                                                                    ))
            except subprocess.CalledProcessError as e:
                string = bashCommand.split(' ')
                print("Genomas: {string1}, Ejecucion k={k} p={p} .-".format(string1=string[1],
                                                                                                                        g=g,
                                                                                                                        k=k,
                                                                                                                        p=p
                                                                                                                    ))