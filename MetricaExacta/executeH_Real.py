import subprocess

#genomas faltantes para ejecutar sobre todos los 17 genomas
genomas = [10, 12, 14, 15, 19, 29, 33, 35, 38, 40, 41, 45, 46, 47, 56, 62, 72, 73, 82, 83]
#topk = [32, 64, 128, 256, 512, 1024, 2048]
topK=128
for g in genomas:
    for k in [10,15]:
        #for value in topk:
        bashCommand = "./entropy_Real GCF_00{g} {k} {value} ".format(g=g, k=k, value=topK)
        try:
            subprocess.run(bashCommand, shell=True, check=True, capture_output=True)
            string = bashCommand.split(' ')
            print("Genomas: {string1}, Ejecucion k={k}, topk={value}.-".format(string1=string[1],
                                                                                  g=g,
                                                                                  k=k,
                                                                                  value=topK))
        except subprocess.CalledProcessError as e:
            string = bashCommand.split(' ')
            print("Genomas: {string1}, Ejecucion k={k}, topk={value}.-".format(string1=string[1],
                                                                                  g=g,
                                                                                  k=k,
                                                                                  value=topK))