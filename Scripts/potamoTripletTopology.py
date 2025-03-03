import dendropy
import sys

def potamoTripletTopology(treefile):
    tree = dendropy.Tree.get(path=treefile,schema='newick',preserve_underscores=True)
    fileSplit = treefile.split('_')
    treeString = tree.as_string(schema='newick')
    pD = {}
    pdc = tree.phylogenetic_distance_matrix()
    for i, t1 in enumerate(tree.taxon_namespace[:-1]):
        for t2 in tree.taxon_namespace[i+1:]:
            pD[(t1.label,t2.label)] = pdc(t1, t2)
            pD[(t2.label,t1.label)] = pdc(t1, t2)
    totalLeaves = []
    for taxon in tree.taxon_namespace:
        totalLeaves.append(taxon.label)
        if 'Pest' in taxon.label:
            peLeaf = taxon.label
        elif 'Pant' in taxon.label:
            paLeaf = taxon.label
        elif 'Pkai' in taxon.label:
            pkLeaf = taxon.label
    tree.reroot_at_midpoint(update_bipartitions=False)
    topology = 'N/A'
    for node in tree:
        subTree = node.leaf_nodes()
        currLeaves = []
        for leaf in subTree:
            currLeaves.append(leaf.taxon.label)
        if len(currLeaves) == 2:
            if paLeaf in currLeaves and pkLeaf in currLeaves:
                topology = 'PaPk'
            elif peLeaf in currLeaves and pkLeaf in currLeaves:
                topology = 'PePk'
            elif peLeaf in currLeaves and paLeaf in currLeaves:
                topology = 'PaPe'
    bootFile = open(fileSplit[0] + '_phyml_boot_trees_' + fileSplit[-1],'r')
    bootLines = bootFile.readlines()
    bootFile.close()
    midpointBoot = 0
    for bootTreeString in bootLines:
        bootTree = dendropy.Tree.get(data=bootTreeString,schema='newick',preserve_underscores=True)
        bootTree.reroot_at_midpoint(update_bipartitions=False)
        bootTopology = False
        for node in bootTree:
            subTree = node.leaf_nodes()
            currLeaves = []
            for leaf in subTree:
                currLeaves.append(leaf.taxon.label)
            if len(currLeaves) == 2:
                if paLeaf in currLeaves and pkLeaf in currLeaves:
                    bootTopology = 'PaPk'
                elif peLeaf in currLeaves and pkLeaf in currLeaves:
                    bootTopology = 'PePk'
                elif peLeaf in currLeaves and paLeaf in currLeaves:
                    bootTopology = 'PaPe'
        if bootTopology == topology:
            midpointBoot += 1
    geneSplit = fileSplit[0].split('.')
    sys.stdout.write(geneSplit[0] + '\t' + geneSplit[3] + '\t' + treeString[0:-1] + '\t' + topology + '\t' + str(midpointBoot) + '\t' + str(pD[(paLeaf,pkLeaf)]) + '\t' + str(pD[(peLeaf,pkLeaf)]) + '\t' + str(pD[(paLeaf,peLeaf)]) + '\n')



sys.stdout.write('Gene\ttrimmingStrategy\tunrootedNewickString\tmidpointRootedTopology\tmidpointBootstrap\tPa_Pk_patristicDistance\tPe_Pk_patristicDistance\tPa_Pe_patristicDistance\n')
fofn = open(sys.argv[1],'r')
for line in fofn:
    potamoTripletTopology(line[0:-1])