import subprocess

genomas = [[83,10], 
           [33,12], 
           [10,19],
           [29,19],
           [40,19]
]
"""            [19,62],
           [29,10],
           [29,40],
           [29,56],
           [35,38],
           [40,10],
           [40,38],
           [41,10],
           [41,19],
           [41,29],
           [41,40],
           [41,79],
           [45,46],
           [47,52],
           [56,62],
           [69,71],
           [72,73],
           [76,79],
           [82,83],
           [89,92],
           [96,98]
""""""  genomas = [[14,15],
           [79,14],
           [35,15],
           [14,69]
           ] """ """ """
pq = [[4, 2], [4, 3], [4, 4], [5, 4], [5, 5], [6, 5]]
for g in range(0, 5):
    for k in [10, 15]:
        g1 = genomas[g][0]
        g2 = genomas[g][1]
        if k == 10:
            for p in range(3,6):
                i = pq[p][0]
                j = pq[p][1]
                bashCommand = "./est_JSD GCF_00{g1}.fna GCF_00{g2}.fna {k} 12 {i} {j} 12 3".format(k=k, g1=g1, g2=g2, i=i, j=j)
                try:
                    subprocess.run(bashCommand, 
                                            shell=True,
                                            check=True, 
                                            capture_output=True)
                except subprocess.CalledProcessError as e:
                    string = bashCommand.split(' ')
                    print("Ejecucion: {string1} {string2}, Ejecucion k={k}, pq_h={i}, pq_w={j} .-".format(string1 = string[1], string2 = string[2], k=k, i=i, j=j))
        else:
            for p in range(0, 3):
                i = pq[p][0]
                j = pq[p][1]
                bashCommand = "./est_JSD GCF_00{g1}.fna GCF_00{g2}.fna {k} 12 {i} {j} 12 3".format(
                    k=k, g1=g1, g2=g2, i=i, j=j)
                try:
                    subprocess.run(bashCommand,
                                shell=True,
                                check=True,
                                capture_output=True)
                except subprocess.CalledProcessError as e:
                    string = bashCommand.split(' ')
                    print("Ejecucion: {string1} {string2}, Ejecucion k={k}, pq_h={i}, pq_w={j} .-".format(
                        string1=string[1], string2=string[2], k=k, i=i, j=j)) 
 