import subprocess

genomas = [[10, 83],
           [12, 33],
           [19, 10],
           [19, 29],
           [19, 40],
           [19, 62],
           [29, 10],
           [29, 40],
           [29, 56],
           [35, 38],
           [40, 10],
           [40, 38],
           [41, 10],
           [41, 19],
           [41, 29],
           [41, 40],
           [45, 46],
           [56, 62],
           [72, 73],
           [82, 83],
           [14, 15],
           [35, 15]
           ]
for k in [10,15]:
    for g in range(0, 22):
        g1 = genomas[g][0]
        g2 = genomas[g][1]
        bashCommand = "./JSD_real GCF_00{g1}.fna GCF_00{g2}.fna {k}".format(g1=g1, g2=g2, k=k)
        try:
            subprocess.run(bashCommand, 
                            shell=True,
                            check=True, 
                            capture_output=True)
            string = bashCommand.split(' ')
            print("Genomas: {string1} | {string2} Ejecucion k={k}".format(string1=string[1], string2=string[2], k=k))
        except subprocess.CalledProcessError as e:
            print("Exception")