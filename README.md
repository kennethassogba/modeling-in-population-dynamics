## Discrétisation d’Equations Différentielles

# Modélisation en dynamique des populations

##Projet
L’objectif de ce travail est de présenter plusieurs modélisations de l’évolution de populations : population isolée, populations composées d’éspèces en “prédation” (prédateurs et de proies), d’éspèces en comptition (partageant la même ressource, ...). Dans chaque cas, la modélisation est donnée par une EDO en dimension 1 ou 2. Pour chaque EDO, vous devez vous poser et répondre aux questions suivantes : quels sont les états d’équilibre ?(s’il y en a). Ces états sont ils stables ? instables ? Selon les cas, vous devez répondre à ces questions en faisant des calculs à la main ou en vous aidant des fonctions de Scilab. Vous devez commenter vos programmes et analyser les résultats obtenus. N’hésitez pas à modifier les paramètres (ainsi que les conditions initiales) donnés dans le texte pour voir et comprendre leur influence sur la solution de l’ODE.

##Introduction
La modélisation des dynamiques des populations vise à expliquer, et éventuellement à prévoir, les évolutions d’une population. La plupart du temps, elle se limite à la description des variations de la taille de la population, mais elle peut aussi permettre de décrire l’évolution de sa structure en âge, voire de sa structure génétique. De nombreux scientifiques ont, proposé des modèles, c’est-à-dire des représentations simplifiées de telles dynamiques,
le plus souvent sous forme mathématique.

Dans ce rapport pour notre projet de discrétisation d’EDO, nous étudierons
quatres modèles notamment ceux de Lotka-Volterra et de Verhulst.

Dans chaque cas, nous ferons une étude des points d’équilibre et une résolution
numérique avec Scilab. Enfin nous présenterons une généralisation de l’un
des modèles.