# Lab notebook - Projet de Physique Numérique 

## 04-11-2026 : Séance 1

### --- Article : Smectic-like bundle formation of planktonic bacteria upon nutrient starvation

Dans cet article (voir le pdf pour plus de détails), les bactéries privées de nutriments se rassemblent en fagots. Cette formation est causée par la sécretion déléments extracellulaire qui, par un effet de dépletion, fait rapporcher les bactéries et les grosse particules. 

### --- But du projet / Objectifs ---

Le but de ce projet est de reconstituer le système et le mettre en équation. Il serait idéal de pouvoir simuler une visualisation d système. 

Parmis les outils à utiliser et maitriser : physique théorique, simulation et résolution d'équation différentielles, visualisation matplotlib, animation, python orienté-objet, calcul dans C/C++ et appréciation de la vitesse de simulation/nombre de paramètre et évolution du système dans sa complexité. 

### --- Discussion avec C. Deroulers ---

On peut séparer le système d'étude en trois différents paramètres à simuler. Chaques paramètres doivent être étudiés et comparer afin de savoir dans quelles conditions on dois résoudre les équations. Parmi elles, le question de la discretisation (ou continuité) dois être posée. En effet, pour des particules discrètes, on aurait un nombre N de particules avec une position $(x,y)$, une vitesse associée. Cependant, pour les systèmes très nombreux mais très petit, on généralise ces particules sous forme d'une densité de particules (ou concentration).  
- Nutriments : $n(x,y,t)$. Il est généralement considéré continu au vu de la taille de ces particules. 
- Extracellular polymerous substances (EPS) : $e(x,y,t)$. Produite pas les cellules, il serait plus compliqué de les considérer continues en concentration car il s'agit de molécules plus grosses que des nutriments mais plus petites que des bactéries.
- Bactéries : $b(x,y,t)$. Elles seront considérée ponctuelles. On ne prendra pas en compte leur multiplication ou la mort cellulaire.

**Pour chaques paramètres, il faut vérifier l'ordre de grandeur**. Il faudra aussi **vérifier les constantes relatives à la diffusion de chaques paramètres**. En effet, si la déplétion de nutriment se fait autours d'une bactérie sur une échelle bien inférieure à la bactérie, alors la concentration en nutriment est la même pour le reste des bactéries. Dans la variation temporelle de la concentration en nutriment : 
$$\frac{\partial n(x,y,t)}{\partial t} = D_n \nabla^2 n(x,y,t) - f_c(x,y,t)$$
Avec $f_c(x,y,t)$ la fonction de "consommation" en nutriment. On devra donc comparer les termes de diffusion et de consommation pour savoir comment évolue les différents paramètres. 

On devra aussi imposer des conditions aux bord, notamment des conditions sur les flux et les concentrations (par exemple, concentration  nulle aux bords). 

**Première étape : simuler à une dimension une source de d'EPS à z=0 pour lesquelles :**

$$\frac{\partial n(x,t)}{\partial t} = D_n \nabla^2 n(x,t)\,,\,n(0, 0)=1$$

#### Un mot sur la méthode d'euler ou le RK4 : 
La méthode d'Euler est une méthode d'intégration numérique qui permet la résolution d'edo. Elle repose sur la décomposition de la dérivé suivant un développement de Taylor : 

$$\frac{dx}{dt}\approx \frac{x(t+\Delta T)-x(t)}{\Delta T}$$

Ici, $dx/dt$ est connu depuis l'edo. On peut faire plusieurs edo consécutives afin d'avoir des edo de degré suppérieurs. On peut développer les développements limité à des ordres suppérieurs : 

$$\frac{d^2x}{dt^2}\approx \frac{x(t+\Delta T)+x(t-\Delta T) - 2x(t)}{\Delta T^2} + O(\Delta T^2)$$

Les méthodes d'Eulers sont rapides mais divergent rapidement pour les EDO successives (voir l'exercice sur la résolution de l'oscillateur). 

La méthode de Runge-Kutta de 4e ordre (RK4) vient apporter plus de stabilité et de précision. Au lieu d'utiliser la dérivé directe dx/dt au temps k, Runge-Kutta prédit la dérivée au centre de l'interval ($y_c(t+\Delta T/2) = y(n) + d_1\Delta T/2$) avec $d_1$ la dérivée au point y(t).  puis on calcule la dérivée au point $y_c$ ($d_2$) et on répète pour avoir $y(t+\Delta T)$. Ça, c'est le RK2, les erreurs sont de l'ordre de $\Delta t^3$ c'est à dire qu'on peut avoir le même résolution qu'un Euleur avec $10^2$ intervals en moins. Cependant, le RK demande plus de calculs donc peut-être plus de temps de computation. 

Remarque : le RK 4 se fait en 4 étapes qui donne une erreur de l'ordre $\Delta t^5$ : 

$k_1 = hf(x_n, y_n)$
$k_2 = hf\left(x_n + \frac{h}{2}, y_n + \frac{k_1}{2}\right)$
$k_3 = hf\left(x_n + \frac{h}{2}, y_n + \frac{k_2}{2}\right)$
$k_4 = hf(x_n + h, y_n + k_3)$

#### Suite de la discusion sur le projet : 

D'après l'article de Takeuchi, les fagots se forment suite à des effets entropiques de depletion entre les bactéries. La formation d'EPS provoque un potentiel de dépletion qui dépends de la position au niveau de la bactérie (extrémité ou coté plat) mais aussi de la concentration en agent de dépletion (EPS). Deux approches peuvent être faites : 
- Approche 1 : interraction singulière entre chaques bactéries du type $U_{1-2}$ entre les bactéries 1 et 2. Il faudra donc sommer l'ensemble des énergies modulo la distance et l'orientation des interactions pour en déduire le vecteru de force et donc le déplacement par un PFD. 
- Approche 2 : Interaction générée comme un champs de potentiels $U$ pour lequel, l'addition de chaques champs mênera à une carte des potentiel décendant et donc aux décentes de gradient nécessaires à toutes les bactéries. 

Ces deux approches font l'hypothèse que les EPS sont des particules bien inférieures en taille que les bactéries et que l'on peut les modélier comme un champs de densité (ou de concentration). 

Dans le cas où modéliser ce champs ne serait pas juste physiquement, on peut toujours les modéliser comme des sphères rigides de rayon $r$ qui se comporterait comme du sable autour de grosse particules, les bactéries. 

### --- Suite des évènements ---

1. Formuler les équations du mouvement et conditions aux bords pour le mouvement de bactéries. Peut être commencé par un mouvement en trajectoire aléatoire et diffusion + collision (type billard avec bruits). 

2. Formuler les équations pour la production d'EPS en un point de l'espace dans l'hypothèse de champs de concentration. Il faudra un point source. 

3. Se renseigner sur la résolution numérique d'EDO avec C/C++ pour accélérer la résolution pour chaques systèmes. 

---
**FIN DE SEANCE**
---
## Interséance 1 

### Création d'un GitHub repository pour stockage de matériel + traçage des différentes versions de codes et du projet

Création du repository sur GitHub avec un README.md et un .gitignore pour ne pas partager les fichier .pdf.

## 09-11-2026 : Séance 2

### --- Premier programme : Résolution d'une équation de diffusion ---

