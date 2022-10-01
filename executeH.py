import subprocess

genomas = [14,10]
pq = [[3, 2],
      [3, 3],
      [4, 2],
      [3, 4],
      [4, 3],
      [4, 4],
      [5, 4],
      [5, 5],
      [6, 4],
      [6, 5]]

for g in genomas:
    if g==10:
        for k in [10,15]:
            if k==15:
                for p in range(0, 10):
                    i = pq[p][0]
                    j = pq[p][1]
                    c1 = [11, 2]
                    c2 = [12, 2]
                    c3 = [12, 3]           
                    bashCommand = "./entropy GCF_00{g} {k} 12 {i} {j} {a} {b}".format(g=g,
                            k=k, i=i, j=j, a=c1[0], b=c1[1])
                    try:
                        subprocess.run(bashCommand, 
                                    shell=True,
                                    check=True, 
                                    capture_output=True)
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
                            
                                
                    bashCommand = "./entropy GCF_00{g} {k} 12 {i} {j} {a} {b}".format(g=g,
                            k=k, i=i, j=j, a=c2[0], b=c2[1])
                    try:
                        subprocess.run(bashCommand, 
                                    shell=True,
                                    check=True, 
                                    capture_output=True)
                        string = bashCommand.split(' ')
                        print("Genomas: {string1}, Ejecucion k={k}, pq_h={i}, pq_w={j}, a={a}, b={b} .-".format(string1=string[1],
                                                                                                                        g=g,
                                                                                                                        k=k,
                                                                                                                        i=i,
                                                                                                                        j=j,
                                                                                                                        a=c2[0],
                                                                                                                        b=c2[1]))
                    except subprocess.CalledProcessError as e:
                        string = bashCommand.split(' ')
                        print("Genomas: {string1}, Ejecucion k={k}, pq_h={i}, pq_w={j}, a={a}, b={b} .-".format(string1=string[1],
                                                                                                                        g=g,
                                                                                                                        k=k,
                                                                                                                        i=i,
                                                                                                                        j=j,
                                                                                                                        a=c2[0],
                                                                                                                        b=c2[1]))
                                
                    bashCommand = "./entropy GCF_00{g} {k} 12 {i} {j} {a} {b}".format(g=g, 
                                                                                            k=k, 
                                                                                            i=i, 
                                                                                            j=j, 
                                                                                            a=c3[0], 
                                                                                            b=c3[1])
                    try:
                        subprocess.run(bashCommand, 
                                    shell=True,
                                    check=True, 
                                    capture_output=True)
                        string = bashCommand.split(' ')
                        print("Genomas: {string1}, Ejecucion k={k}, pq_h={i}, pq_w={j}, a={a}, b={b} .-".format(string1=string[1],
                                                                                                                        g=g,
                                                                                                                        k=k,
                                                                                                                        i=i,
                                                                                                                        j=j,
                                                                                                                        a=c3[0],
                                                                                                                        b=c3[1]))

                    except subprocess.CalledProcessError as e:
                        string = bashCommand.split(' ')
                        print("Genomas: {string1}, Ejecucion k={k}, pq_h={i}, pq_w={j}, a={a}, b={b} .-".format(string1=string[1],
                                                                                                                        g=g,
                                                                                                                    k=k,
                                                                                                                    i=i,
                                                                                                                    j=j,
                                                                                                                    a=c3[0],
                                                                                                                    b=c3[1]))
