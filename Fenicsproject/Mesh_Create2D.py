import dolfin as df


def create_mesh(txt):
    # pegando pontos do arquivo vtu e definindo nossos numeros de vertices (nvertices) e o numero de celulas (ncells)
    with open(txt) as arquivo:
        ret = arquivo.read()
        ret = ret.split('"')
    nvertices = int(ret[7])
    ncells = int(ret[9])
    # limpando o arquivo.vtu de caracteres especiais e deixando ele com todos os valor como elementos da lista limpa
    with open(txt) as arquivo:
        ret = arquivo.read()
        ret = ret.split("'", )

    nova = []
    novamente = []
    limpa = []
    for palavras in ret:
        nova.append(palavras.split('\n'))
    for i in range(len(nova)):
        for palavras in nova[i]:
            novamente.append(palavras.split(' '))

    for i in range(len(novamente)):
        for palavras in novamente[i]:
            if palavras is not '':
                limpa.append(palavras)

    print(f'o numero de vertices é {nvertices} \n')
    print(f'o numero de faces é {ncells} \n')
    # criando as listas de vertices e de celulas, a lista de vertices começa no indice 13 pois é o indice que começa os pontos dela, indo até 3 vezes o numero de vertices,pois cada vertice tem x,y e z, as celular começam 7 indices depois do fim do numero de vertices indo até 3 vezes o numero dela pois tambem tem 3 elementos pra cada vertice
    vertices = [float(limpa[i]) for i in range(13, 13 + 3 * nvertices)]
    cells = [int(limpa[i]) for i in range(20 + 3 * nvertices, 20 + 3 * nvertices + 3 * ncells)]

    # para poder criar a mesh foi criado um vetor de pontos dolfin utilizando a lista vertices (Vpoints), cada coordenadas de vertices é criado um ponto dolfin, já para as celulas era necessario um array com 3 elementos, então foi criado uma lista onde cada elemento da lista é um array contendo 3 coordenadas da celula(Cpoints)
    Vpoints = []
    for i in range(0, len(vertices), 3):
        Vpoints.append(df.Point(vertices[i], vertices[i + 1]))
    Cpoints = []
    for i in range(0, len(cells), 3):
        Caux = [cells[i], cells[i + 1], cells[i + 2]]
        Cpoints.append(Caux)

    # criando a malha
    mesh = df.Mesh()
    editor = df.MeshEditor()
    editor.open(mesh, 'triangle', 2, 2)

    # Adicionandos os vertices contidos em Vpoints
    editor.init_vertices(nvertices)
    for i in range(nvertices):
        editor.add_vertex(i, Vpoints[i])

    # adicionando as celulas que estão em Cpoints
    editor.init_cells(ncells)
    for i in range(len(Cpoints)):
        editor.add_cell(i, Cpoints[i])
    # fechando a edição da malha e retornando a mesma
    editor.close()
    print(f'criado uma malha com {nvertices} vertices e {ncells} celulas')

    return mesh

