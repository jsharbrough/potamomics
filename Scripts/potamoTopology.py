import dendropy
import sys

def potamoTopology(treefile):
	try:
		tree = dendropy.Tree.get(path=treefile,schema='newick',preserve_underscores=True)
	except:
		sys.stderr.write(treefile + ' not found\n')
		return
	fileSplit = treefile.split('.')
	og = fileSplit[-3]
	pD = {}
	pdc = tree.phylogenetic_distance_matrix()
	for i, t1 in enumerate(tree.taxon_namespace[:-1]):
		for t2 in tree.taxon_namespace[i+1:]:
			pD[(t1.label,t2.label)] = pdc(t1, t2)
	 		pD[(t2.label,t1.label)] = pdc(t1, t2)
	totalLeaves = []
	paLeaves = []
	for taxon in tree.taxon_namespace:
		totalLeaves.append(taxon.label)
		if 'Pest' in taxon.label:
			peLeaf = taxon.label
		elif 'Pant' in taxon.label:
			paLeaves += [taxon.label]
		elif 'Pkai' in taxon.label:
			pkLeaf = taxon.label
	sys.stdout.write(og + '\t' + str(pD[(paLeaves[0],paLeaves[1])]) + '\t' + str(pD[(paLeaves[0],pkLeaf)]) + '\t' + str(pD[(paLeaves[1],pkLeaf)]) + '\t' + str((pD[(paLeaves[0],pkLeaf)] + pD[(paLeaves[1],pkLeaf)])/2.0) + '\t' + str(pD[(paLeaves[0],peLeaf)]) + '\t' + str(pD[(paLeaves[1],peLeaf)]) + '\t' + str((pD[(paLeaves[0],peLeaf)] + pD[(paLeaves[1],peLeaf)])/2.0) + '\t' + str(pD[(pkLeaf,peLeaf)]) + '\n')
	'''tree.reroot_at_midpoint(update_bipartitions=False)
	treeString = tree.as_string(schema='newick')
	topology = 'N/A'
	bootValue = 'N/A'
	for node in tree:
		subTree = node.leaf_nodes()
		currLeaves = []
		for leaf in subTree:
			currLeaves.append(leaf.taxon.label)
		if len(currLeaves) == 2:
			if paLeaf in currLeaves and pkLeaf in currLeaves:
				topology = 'PaPk'
				bootValue = node.label
			elif peLeaf in currLeaves and pkLeaf in currLeaves:
				topology = 'PePk'
				bootValue = node.label
			elif peLeaf in currLeaves and paLeaf in currLeaves:
				topology = 'PaPe'
				bootValue = node.label
	sys.stdout.write(fileSplit[0] + '\t' + treeString[0:-1] + '\t' + topology + '\t' + str(bootValue) + '\t' + str(pD[(paLeaf,pkLeaf)]) + '\t' + str(pD[(peLeaf,pkLeaf)]) + '\t' + str(pD[(paLeaf,peLeaf)]) + '\n')'''
	


sys.stdout.write('og\tPaPa_PatristicDistance\tPa1Pk_PatristicDistance\tPa2Pk_PatristicDistance\tmeanPaPk_PatristicDistance\tPa1Pe_PatristicDistance\tPa2Pe_PatristicDistance\tmeanPaPe_PatristicDistance\tPkPe_PatristicDistance\n')
fofn = open(sys.argv[1],'r')
for line in fofn:
	realLine = line
	while realLine[-1] == '\t' or realLine[-1] == '\r' or realLine[-1] == '\n':
		realLine = realLine[0:-1]
	potamoTopology(realLine)