import subprocess

#genomas faltantes para ejecutar sobre todos los 17 genomas
#genomas = [12, 19, 29, 33, 35, 38, 40, 41, 45, 46, 47, 56, 62, 72, 73, 82, 83]
genomas = [14, 15]
pq = [[4,2], [4,3], [4, 4], [5, 4], [5, 5], [6, 5]]
#k=10
for g in genomas:
    for k in [10,15]:
        if k==10:
            for p in range(3, 6):
                    i = pq[p][0]
                    j = pq[p][1]
                    c1 = [12, 3]
                    bashCommand = "./entropy GCF_00{g} {k} 12 {i} {j} {a} {b}".format(g=g, k=k, i=i, j=j, a=c1[0], b=c1[1])
                    try:
                        subprocess.run(bashCommand, shell=True, check=True, capture_output=True)
                        string = bashCommand.split(' ')
                        print("Genomas: {string1}, Ejecucion k={k}, pq_h={i}, pq_w={j}, a={a}, b={b} .-".format(string1=string[1],
                                                                                                                            g=g,
                                                                                                                            k=k,
                                                                                                                            i=i,
                                                                                                                            j=j,
                                                                                                                            a=c1[0],
                                                                                                                            b=c1[1]))
                    except subprocess.CalledProcessError as e:
                        string = bashCommand.split(' ')
                        print("Genomas: {string1}, Ejecucion k={k}, pq_h={i}, pq_w={j}, a={a}, b={b} .-".format(string1=string[1],
                                                                                                                            g=g,
                                                                                                                            k=k,
                                                                                                                            i=i,
                                                                                                                            j=j,
                                                                                                                            a=c1[0],
                                                                                                                            b=c1[1]))
        else:
            for p in range(0, 3):
                    i = pq[p][0]
                    j = pq[p][1]
                    c1 = [12, 3]
                    bashCommand = "./entropy GCF_00{g} {k} 12 {i} {j} {a} {b}".format(g=g, k=k, i=i, j=j, a=c1[0], b=c1[1])
                    try:
                        subprocess.run(bashCommand, shell=True, check=True, capture_output=True)
                        string = bashCommand.split(' ')
                        print("Genomas: {string1}, Ejecucion k={k}, pq_h={i}, pq_w={j}, a={a}, b={b} .-".format(string1=string[1],
                                                                                                                            g=g,
                                                                                                                            k=k,
                                                                                                                            i=i,
                                                                                                                            j=j,
                                                                                                                            a=c1[0],
                                                                                                                            b=c1[1]))
                    except subprocess.CalledProcessError as e:
                        string = bashCommand.split(' ')
                        print("Genomas: {string1}, Ejecucion k={k}, pq_h={i}, pq_w={j}, a={a}, b={b} .-".format(string1=string[1],
                                                                                                                            g=g,
                                                                                                                            k=k,
                                                                                                                            i=i,
                                                                                                                            j=j,
                                                                                                                            a=c1[0],
                                                                                                                            b=c1[1]))