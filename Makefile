FLAGS=-std=c++17 -Wall -Werror 

#dossier de compilation pour centraliser les fichiers en .o
BUILD_DIR=./build

#liste des fichiers de chaque version (fichier main exclus)
BUILD_PRELIM=research Protein Database Index Utils Tas_min
BUILD_FINAL=$(BUILD_PRELIM)
BUILD_OPT= researchopt Protein Database Index Utils Tas_min research

#objets de chaque version à créer dans le dossier build
OBJ_PRELIM=$(patsubst %,$(BUILD_DIR)/%.o,$(BUILD_PRELIM))
OBJ_FINAL=$(patsubst %,$(BUILD_DIR)/%.o,$(BUILD_FINAL))
OBJ_OPT=$(patsubst %,$(BUILD_DIR)/%.o,$(BUILD_OPT))

#recette de création des fichiers .o dans le dossier build
$(BUILD_DIR)/%.o: %.cpp
	echo $@
	g++ -c $^ -o $@ $(FLAGS)

#version préliminaire
.PHONY: projetprelim
projetprelim: build_dir projetprelim.cpp $(OBJ_PRELIM)
	g++ -o $@ $@.cpp $(OBJ_PRELIM) ${FLAGS} 

#version finale
.PHONY: projet
projet: build_dir projet.cpp $(OBJ_FINAL)
	g++ -o $@ $@.cpp $(OBJ_FINAL) ${FLAGS}

#version optimisée
.PHONY: projetopt
projetopt: build_dir projetopt.cpp $(OBJ_OPT)
	g++ -o $@ $@.cpp $(OBJ_OPT) ${FLAGS}

#crée un dossier build qui stocke les fichiers .o
build_dir:
	mkdir -p $(BUILD_DIR)

#supprime le dossier build et les fichiers dedans
clean: 
	rm -rf $(BUILD_DIR)
