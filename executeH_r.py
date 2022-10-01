import subprocess

genomas = [10, 14, 15]

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
    for k in [10, 15]:
        for p in range(0, 10):
            i = pq[p][0]
            j = pq[p][1]
            bashCommand = "./entropy_CR GCF_00{g} {k} 12 {i} {j}".format(
                g=g, k=k, i=i, j=j)
            try:
                subprocess.run(bashCommand, shell=True,
                               check=True, capture_output=True)
                string = bashCommand.split(' ')
                print("Genoma: {string1}, Ejecucion k={k}, pq_h={i}, pq_w={j}.-".format(string1=string[1],
                                                                                                        g=g,
                                                                                                        k=k,
                                                                                                        i=i,
                                                                                                        j=j))
            except subprocess.CalledProcessError as e:
                string = bashCommand.split(' ')
                print("Genoma: {string1}, Ejecucion k={k}, pq_h={i}, pq_w={j}.-".format(string1=string[1],
                                                                                                        g=g,
                                                                                                        k=k,
                                                                                                        i=i,
                                                                                                        j=j))