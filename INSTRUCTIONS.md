# Instructions de Compilation et d'Exécution

Ce projet contient des solveurs pour l'équation différentielle `-u'' = f` (en 1D) et `-Δu = f` (en 2D) en utilisant les méthodes des Différences Finies et des Volumes Finis.

## Structure du Projet

Le code source Java se trouve dans `src/main/java/`. Chaque méthode a son propre package :
*   `FiniteDifference1` : Différences Finies 1D
*   `FiniteVolume1` : Volumes Finis 1D
*   `FiniteDifference2` : Différences Finies 2D
*   `FiniteVolume2` : Volumes Finis 2D
*   `LaplaceFiniteDifference2D` : Solveur pour l'équation de Laplace en 2D (`-Δu = 0`) avec `u(x,y)=2x+y` sur maillage variable.

Chaque package contient une classe principale avec une méthode `main` qui lance une interface graphique Swing.

## Prérequis

*   JDK (Java Development Kit) installé (par exemple, version 8 ou ultérieure).
*   Accès à un terminal ou une ligne de commande.

## Compilation

1.  Ouvrez un terminal ou une ligne de commande.
2.  Naviguez jusqu'à la racine du projet (où se trouve le dossier `src` et le `pom.xml`).
3.  Compilez tous les fichiers Java. Étant donné la présence d'un `pom.xml`, Maven peut être utilisé, mais une compilation manuelle est aussi simple pour ce projet sans dépendances externes complexes autres que Swing (inclus dans le JDK).

    **Option A : Compilation manuelle (si pas de Maven ou pour simplicité)**
    Créez un répertoire pour les classes compilées (par exemple, `out`):
    ```bash
    mkdir out
    ```
    Compilez tous les fichiers `.java` en spécifiant le répertoire de sortie et l'encodage UTF-8 pour une gestion correcte des accents. Depuis la racine du projet :
    ```bash
    javac -encoding UTF-8 -d out src/main/java/FiniteDifference1/ODEFiniteDifference.java src/main/java/FiniteVolume1/ODEFiniteVolume.java src/main/java/FiniteDifference2/ODEFiniteDifference.java src/main/java/FiniteVolume2/ODEFiniteVolume.java src/main/java/LaplaceFiniteDifference2D/LaplaceSolverGUI.java
    ```
    Si vous utilisez Java 11+ et que les classes `HeatmapPanel` sont laissées comme classes de premier niveau (non-statiques imbriquées ou publiques séparées) dans les fichiers 2D, vous pourriez avoir besoin de les compiler explicitement ou de les rendre statiques imbriquées / publiques. Le code fourni les a comme classes de premier niveau dans le même fichier, ce qui est valide.

    **Option B : Utilisation de Maven (recommandé si configuré)**
    Si Maven est installé et que le `pom.xml` est configuré pour la compilation (même simple), vous pouvez généralement compiler avec :
    ```bash
    mvn compile
    ```
    Les classes compilées se trouveront alors dans le répertoire `target/classes/`. (Note: le `pom.xml` actuel ne semble pas configuré pour un build Java standard, il faudra peut-être l'adapter ou utiliser la compilation manuelle).

## Exécution

Après la compilation, vous pouvez exécuter chaque programme en spécifiant le classpath et la classe principale.

**Si compilé manuellement dans `out` (depuis la racine du projet) :**

*   **Différences Finies 1D:**
    ```bash
    java -cp out FiniteDifference1.ODEFiniteDifference
    ```
*   **Volumes Finis 1D:**
    ```bash
    java -cp out FiniteVolume1.ODEFiniteVolume
    ```
*   **Différences Finies 2D:**
    ```bash
    java -cp out FiniteDifference2.ODEFiniteDifference
    ```
*   **Volumes Finis 2D:**
    ```bash
    java -cp out FiniteVolume2.ODEFiniteVolume
    ```
*   **Solveur Laplace 2D avec UI (Maillage Variable):**
    ```bash
    java -cp out LaplaceFiniteDifference2D.LaplaceSolverGUI
    ```

**Si compilé avec Maven (depuis la racine du projet, après `mvn compile`) :**
Le classpath est `target/classes/`.

*   **Différences Finies 1D:**
    ```bash
    java -cp target/classes FiniteDifference1.ODEFiniteDifference
    ```
*   **Volumes Finis 1D:**
    ```bash
    java -cp target/classes FiniteVolume1.ODEFiniteVolume
    ```
*   **Différences Finies 2D:**
    ```bash
    java -cp target/classes FiniteDifference2.ODEFiniteDifference
    ```
*   **Volumes Finis 2D:**
    ```bash
    java -cp target/classes FiniteVolume2.ODEFiniteVolume
    ```

## Utilisation des Interfaces Graphiques

*   **Programmes 1D (`ODEFiniteDifference` dans `FiniteDifference1`, `ODEFiniteVolume` dans `FiniteVolume1`) :**
    *   Une fenêtre s'ouvrira.
    *   **Sélectionnez la "Solution exacte u(x)" de référence** (pour laquelle `f(x) = -u''(x)` sera calculé). Les options sont `u(x) = sin(πx)` et `u(x) = x³`.
    *   **Sélectionnez le "Type d'équation"**. Pour cette tâche, utilisez toujours `TYPE3 ("-u'' = f")`.
    *   Les conditions aux limites sont fixées dans le code à `u(0)=0` et `u(1)=1`.
    *   Cliquez sur **"Résoudre et Afficher Solution"** pour voir des graphiques de la solution numérique vs. analytique pour N=10, 20, 40, 80. Chaque graphique s'ouvre dans un onglet. La console affichera l'erreur L∞.
    *   Cliquez sur **"Afficher Courbe d'Erreur"** pour calculer les solutions pour N=10, 20, 40, 80, 160, 320. La console affichera l'erreur L∞ et l'ordre de convergence. Une nouvelle fenêtre s'ouvrira avec un graphique log-log de l'erreur L∞ en fonction de N.

*   **Programmes 2D (`ODEFiniteDifference` dans `FiniteDifference2`, `ODEFiniteVolume` dans `FiniteVolume2`) :**
    *   Une fenêtre s'ouvrira.
    *   **Sélectionnez la "Sol. exacte u(x,y)" de référence** (pour laquelle `f(x,y) = -Δu(x,y)` sera calculé). Options : `u(x,y) = sin(πx)sin(πy)` et `u(x,y) = x³y³`.
    *   Les conditions aux limites de Dirichlet sont prises à partir de la solution exacte choisie sur les bords du domaine `(0,1)x(0,1)`.
    *   Cliquez sur **"Résoudre et Afficher"**.
    *   La zone de texte affichera l'erreur L∞ et l'ordre de convergence pour N x N avec N = 10, 20, 40, 80.
    *   Un sélecteur "Afficher N:" permet de choisir la taille de maillage pour laquelle les heatmaps sont affichées.
    *   La partie inférieure affichera trois heatmaps pour le N sélectionné : solution numérique, analytique, et erreur absolue.

*   **Solveur Laplace 2D (`LaplaceSolverGUI` dans `LaplaceFiniteDifference2D`) :**
    *   Résout `-Δu = 0` avec `u(x,y) = 2x+y` comme conditions aux limites et solution exacte.
    *   Entrez le "Nombre d'intervalles N" désiré (par exemple, 10, 20, etc.).
    *   Cliquez sur "Résoudre".
    *   La zone de texte affiche l'erreur L∞ pour le N choisi. Si N <= 5, la grille numérique est aussi affichée en texte.
    *   Trois heatmaps (Numérique, Analytique, Erreur Absolue) sont affichées pour le N choisi.

## Interprétation des Résultats

*   **Erreur L∞ :** Différence maximale absolue entre solution numérique et analytique. Plus petite = mieux.
*   **Ordre de Convergence Numérique :** Théoriquement 2 pour les schémas utilisés. Si N double, l'erreur L∞ devrait être divisée par 4.
*   **Graphiques (1D) :**
    *   Solution : Numérique (rouge) vs Analytique (bleu). Devraient se superposer pour N grand.
    *   Courbe d'erreur (log-log) : Droite de pente ~-2 si ordre 2.
*   **Heatmaps (2D) :**
    *   Numérique et Analytique : Devraient se ressembler. Couleurs : Bleu (min) -> Rouge (max).
    *   Erreur Absolue : Devrait être majoritairement bleue (valeurs faibles).

## Note sur le `pom.xml`
Le `pom.xml` présent à la racine semble être un placeholder ou configuré pour un autre type de projet (peut-être lié à `jenny.exe`). Pour compiler les solveurs d'EDP Java, la compilation manuelle avec `javac` comme décrit ci-dessus est la méthode la plus directe si le `pom.xml` n'est pas adapté. Si vous souhaitez utiliser Maven, le `pom.xml` nécessiterait des ajustements pour un projet Java standard (définition des `sourceDirectory`, `outputDirectory`, plugin de compilation Java, etc.).
