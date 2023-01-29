# Chemin relatif

projects/VS2017/raylib.sln

# Remarques

## Limitation

Nous avons rencontré des difficultés à retranscrire certaines équations en c++.
Nous avons eu du mal à obtenir un effet de rotation convaincant.

## Bugs

Rarement, la sphère a des comportements inattendus.

## Voies d'amélioration

Implémenter plus de physique afin de rendre le moteur plus réaliste, par exemple :
- Prendre en compte les frottements d'Euler (de l'air)
- Prendre en compte le mouvement de l'obstacle lors d'une collision
- Optimisation du rendu en affichant que les primitives visibles par la caméra

Ajouter une interface permettant à l'utilisateur de modifier les paramètres de la physique et d'ajouter facilement des primitives 3D, afin de créer l'environnement de test qu'il souhaite.

# Répartiton des tâches

## TD 1 : Primitives 3D & Méthodes de traçage OpenGL associées (polygons & wireframes)

### Caméra

Maxime, Alex et Erwan.

### Conversion

Nous avons fait les conversions à trois.

### Draw

|                 	| Alex 	| Maxime 	| Erwan 	|
|-----------------	|------	|--------	|-------	|
| Box             	|   X  	|    X   	|       	|
| Disk            	|      	|        	|   X   	|
| Sphere          	|   X  	|        	|       	|
| Sphere corner   	|   X  	|        	|       	|
| Hemishpere      	|   X  	|        	|       	|
| Cylindre        	|   X  	|        	|       	|
| Cylinder corner 	|   X  	|        	|       	|
| Capsule         	|   X  	|        	|       	|
| RoundedBox      	|      	|        	|   X   	|
| Segment         	|   X  	|    X   	|       	|

## TD 2 : Méthodes de calcul d’intersection entre un segment et divers objets 3D

### Intersections

|                   	| Alex 	| Maxime 	| Erwan 	|
|-------------------	|------	|--------	|-------	|
| Segment-plane     	|      	|        	|   X   	|
| Quad              	|   X  	|        	|   X   	|
| Disk              	|      	|        	|   X   	|
| Sphere            	|   X  	|        	|       	|
| Box               	|   X  	|        	|       	|
| Infinite cylinder 	|      	|        	|   X   	|
| Cylinder          	|   X  	|        	|   X   	|
| Capsule           	|      	|        	|   X   	|
| RoundedBox        	|   X  	|        	|       	|

## TD 3 : Détection de collision d’une Sphere en mouvement et d’un ensemble de Box / RoundedBox statiques

### Collisions

Vues avec l'ensemble du groupe.

## Interface graphique

Alex s'est bien amusé.
