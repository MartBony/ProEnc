# Tuto Simplifié
## Installation pas à pas
⚠️ Le fichier PDF peut contenir des caractères "entrée" au niveau des cassures de lignes, qui peuvent être interprétées par le terminal comme l'execution d'une commande

($\to$ supprimer les fichiers de l'ancienne configuration de R sur bird si configurée précédemment)

1) Installer MobaXTerm
Puis se connecter à Bird avec son mdp (Session>SSH)
URL : bird2login.univ-nantes.fr
MDP : tonmotdepassebird

(\+ Configurer le raccourci clavier coller dans les paramètres car pratique)

2) Initialiser serveur
```Bash
qlogin
```

3) Installer conda
```Bash
cd /CONDAS/users/mgregorio # A remplacer par id
# Telecharger miniforge (conda)
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
# Donner les permission d'execution
chmod 777 ./Miniforge3-Linux-x86_64.sh
# Executer l'installateur

# ⚠️ spécifier le Répertoire d'installation : /CONDAS/users/monnomdutilisateur/miniforge3
./Miniforge3-Linux-x86_64.sh # yes + Changer le répertoire d'installation + no
```

4) Initialiser conda et ajouter conda au PATH
```bash
export PATH="/CONDAS/users/monnomdutilisateur/miniforge3/bin:$PATH" #Modification temporaire du path pour ajouter la commande conda au système
conda init # Initialisation de conda

cd ~ # Dans le dossier /home/monnomdutilisateur/
nano .bashrc
```

$\to$ Vérifier que le fichier contient un export path avec miniforge3, sinon ajouter le code suivant a la fin du fichier.
⚠️ Ne pas oublier de modifier le nom d'utilisateur
```bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/CONDAS/users/monnomdutilisateur/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/CONDAS/users/monnomdutilisateur/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/CONDAS/users/monnomdutilisateur/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/CONDAS/users/monnomdutilisateur/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/CONDAS/users/monnomdutilisateur/miniforge3/etc/profile.d/mamba.sh" ]; then
    . "/CONDAS/users/monnomdutilisateur/miniforge3/etc/profile.d/mamba.sh"
fi

# <<< conda initialize <<<
```
> Info : .bashrc est exécuté à l'ouverture du terminal et permet de rajouter conda au PATH automatiquement

**Puis fermer l'onglet et rouvrir une session BIRD** + (qlogin)

5) Initialiser environnement
```Bash
cd /CONDAS/users/mgregorio
nano seurenv.yml
```
Puis copier coller dans nano le contenu du fichier seurenv.yml sur le repo.
Sauvegarder (ctrl X + yes)

$\to$ Sortir et enregistrer le fichier

Initialiser l'environnement :
```bash
mamba env create -f seurenv.yml # mamba plus rapide que conda
```
⚠️ C'est long

Puis
```bash
conda env list # Vérifier que seurat_env est installé
```

## Routine d'ouverture
1) Activer l'environnement conda
```bash
conda activate seurat_env
```
2) Ouvrir rstudio
Pour travailler plus proprement, je met mes codes dans un dossier
```bash
mkdir Rcode
```
Ouvrir Rstudio
```bash
rstudio
# Ou avec le nom d'un fichier
rstudio Rcode/monfichier.R
```
$\to$ ignorer les mises à jour à l'ouverture

3) Fermer
Désactiver l'environnement conda
```bash
conda deactivate
```

## Commandes additionnelles
### Mettre à jour l'environement
Après modification du fichier seurenv.yml
```bash
mamba env update --file seurenv.yml --prune
```
> Exemple : pour ajouter un package

### Se déconnecter de l'environement
```bash
conda deactivate
```
### Supprimer l'environnement
Après s'être déconnecté de l'environnement,
```bash
conda remove --name seurat_env --all
```

## Importer ses fichiers sur bird
Utiliser l'interface graphique de mobaXterm par glisser déposer. (Icone planète)