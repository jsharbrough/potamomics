import dendropy
import sys

def potamoQuartetTopology(treefile):
	tree = dendropy.Tree.get(path=treefile,schema='newick',preserve_underscores=True)
	fileSplit = treefile.split('.')
	treeString = tree.as_string(schema='newick')
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
			paLeaves.append(taxon.label)
		elif 'Pkai' in taxon.label:
			pkLeaf = taxon.label
	tree.is_rooted = False
	treeFilter = False
	failingReasons = []
	homoeologs = 'N/A'
	antipodarumBoot = 'N/A'
	for node in tree:
		subTree = node.leaf_nodes()
		currLeaves = []
		for leaf in subTree:
			currLeaves.append(leaf.taxon.label)
		if len(currLeaves) == 2:
			if paLeaves[0] in currLeaves and paLeaves[1] in currLeaves:
				homoeologs = 'homoeologs group'
				if node.label == None or antipodarumBoot == 'N/A':
					antipodarumBoot = node.label
			elif peLeaf in currLeaves and pkLeaf in currLeaves:
				homoeologs = 'homoeologs group'
				if antipodarumBoot == 'N/A' or antipodarumBoot == None:
					antipodarumBoot = node.label
			elif pkLeaf in currLeaves and paLeaves[0] in currLeaves:
				treeFilter = True
				homoeologs = 'homoeologs do not group'
			elif pkLeaf in currLeaves and paLeaves[1] in currLeaves:
				treeFilter = True
				homoeologs = 'homoeologs do not group'
			elif peLeaf in currLeaves and paLeaves[0] in currLeaves:
				treeFilter = True
				homoeologs = 'homoeologs do not group'
			elif peLeaf in currLeaves and paLeaves[1] in currLeaves:
				treeFilter = True
				homoeologs = 'homoeologs do not group'
	topology = 'N/A'
	rootPlacement = 'N/A'
	pa1_pk = pD[(paLeaves[0],pkLeaf)]
	pa2_pk = pD[(paLeaves[1],pkLeaf)]
	pa1_pe = pD[(paLeaves[0],peLeaf)]
	pa2_pe =  pD[(paLeaves[1],peLeaf)]
	pk_pe =  pD[(pkLeaf,peLeaf)]
	pa_pa = pD[(paLeaves[0],paLeaves[1])]
	meanPaPk = (pa1_pk + pa2_pk)/2
	meanPaPe = (pa1_pe + pa2_pe)/2
	if treeFilter == False:
		tree.reroot_at_midpoint(update_bipartitions=False)
		asymmetric = False
		for node in tree:
			subTree = node.leaf_nodes()
			currLeaves = []
			for leaf in subTree:
				currLeaves.append(leaf.taxon.label)
			if len(currLeaves) == 3:
				asymmetric = True
				if paLeaves[0] in currLeaves and paLeaves[1] in currLeaves:
					if pkLeaf in currLeaves:
						topology = 'PaPk'
						rootPlacement = 'Pe'
					elif peLeaf in currLeaves:
						topology = 'PaPe'
						rootPlacement = 'Pk'
				else:
					topology = 'PkPe'
					rootPlacement = 'Pa_singleHomoeolog'
		if asymmetric == False:
			topology = 'PkPe'
			rootPlacement = 'middle'
	og = fileSplit[1]
	sys.stdout.write(og + '\t' + treeString[0:-1] + '\t' + homoeologs + '\t' + str(antipodarumBoot) + '\t' + rootPlacement + '\t' + topology + '\t' + str(round(pa_pa,4)) + '\t' + str(round(meanPaPk,4)) + '\t' + str(round(pk_pe,4)) + '\t' + str(round(meanPaPe,4)) + '\n')
	


sys.stdout.write('Orthogroup_ID\tunrootedNewickString\tHomoeologous_Relationship\tanitpodarumBootstrapSupport\trootPlacement\ttopology\tPa_Pa_patristicDistance\tmeanPa_Pk_patristicDistance\tPk_Pe_patristicDistance\tmeanPa_Pe_patristicDistance\n')
fofn = open(sys.argv[1],'r')
for line in fofn:
	potamoQuartetTopology(line[0:-1])