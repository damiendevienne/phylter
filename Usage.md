  R CMD INSTALL .

  Ensuite, depuis ce dépôt, le CLI s’utilise directement, sans taper Rscript :

  ./exec/phylter --help

  Le jeu de test Carnivora est inclus sous forme de données R ; exporte-le une fois en Newick :

  Rscript -e 'data(carnivora, package="phylter"); ape::write.tree(carnivora, file="carnivora-test.nwk", digits=17)'

  Puis lance l’analyse :

  ./exec/phylter --trees carnivora-test.nwk --out carnivora-result --report

  Tu obtiendras notamment :

  carnivora-result.outliers.tsv
  carnivora-result.discarded.tsv
  carnivora-result.summary.txt
  carnivora-result.session.txt
  carnivora-result.pdf

  Pour un jeu réel avec un fichier Newick multi-arbres :

  ./exec/phylter --trees mes_arbres.nwk --out mon_analyse

  Ou avec un répertoire contenant un arbre Newick par gène :

  ./exec/phylter --trees arbres/ --out mon_analyse --report

  À ce stade, le CLI est plus agréable à appeler, mais il utilise encore R en interne : R CMD INSTALL . et les dépendances R restent nécessaires. L’objectif d’un unique binaire
  autonome, installable comme phyml, demande la réécriture du cœur en C++.
