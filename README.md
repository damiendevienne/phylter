#Phylter is a tool for analyzing, visualizing and filtering Phylogenomics datasets. 

SimOutliersHGT(tree, nbgn, outgn, outsp, outcell, sp = 0)
SimOutliersLg(tree, nbgn, outgn, outsp, outcell, sp = 0)

tree = Arbre d'espèce
nbgn = nombre de gènes souhaité
outgn = nombre d'outliers gènes souhaité
outsp = nombre d'outliers espèces souhaité
outcell = nombre d'outliers cell souhaité
sp = 1 si HGT sur les branches externes uniquement / 0 sinon

	-VIZ : Code d'une méthode de visualisation des gènes et des espèces codé par Damien.

	-PHYLTER : prémisses du package PhylteR. Le code est dans PhylteR/R

	-CLUSTER :

		- ARBRES : Ensemble des listes arbres (numéroté de 1 à 5550). Ce sont des arbres de 100 gènes à 30 espèces. Provenant de 3 arbres d'espèces différents.
			- plan.csv = a quoi correspondent les liste d'arbres en fonction de leurs numéro. Arbres d'espèces => tree1 = ArbreSym.phy / tree2 = ArbreAsym.phy / tree3 = ArbreRandom.phy
			- Les dossiers tree0 à tree4 contiennen les listes d'arbres de 1 à 4500. C'est a dire les arbres contenant ou non des outliers complets.
			- Le dossier treeOutcell ne contient que les listes d'arbres numérotées de 4501 à 5555 contenant les outliers cell (et pas d'outliers complets)

		-SCRIPT : Différents scripts ayant servis à générer les listes d'arbres simulées sur le cluster
		
		-FORETS : jeux de données réèls
		

--> Pour éditer le code du shiny, penser à demander l'acces au svn à Aurélie ou à Stéphane.

