# MemoriaCristopherOsorio


# Neural Branching for Ibex Optimization

Este repositorio contiene la implementación de una arquitectura híbrida para la optimización global de intervalos. Se realiza una modificación en la librería **Ibex** para actuar como cliente TCP y consultar a una **Red Neuronal (MLP)** (Python) que actúa como servidor, prediciendo la mejor heurística de bisección en tiempo de ejecución.

```bash
sudo apt-get update
sudo apt-get install build-essential git bison flex liblapack-dev libblas-dev


git clone [https://github.com/TU_USUARIO/TU_REPOSITORIO.git](https://github.com/TU_USUARIO/TU_REPOSITORIO.git)
cd TU_REPOSITORIO
```

Se debe realizar un clone de la libreria https://github.com/ibex-team/ibex-lib .

En el repoditorio se utilizará el archivo ''ibex_optimizer_conexion.cpp''
El contenido de este archivo será cambiado por el de ''ibex_optimizer.cpp''(recomiendo salvar el contenido del archivo original en otro archivo)



Una vez teniendo el entorno listo, siga las instrucciones del repositorio original para realizar las instalaciones https://ibex-team.github.io/ibex-lib/install-cmake.html :

Luego de intalar, reemplazar de la misma manera el archivo ''ibexopt'' por el de este repositorio, salvando el contenido del archivo original

Luego, ejecutar el archivo ''comunicacionSockets.py'' en comunicación:
```bash
python3 comunicacionSockets.py

```
Luego, ejecutar el optimizer de la siguiente manera:

```bash
$ cd bin
$ mv ibexopt ibex_sockets
$ cd ..
$ bin/ibex_sockets ../benchs/optim/benchs_a_probar
```