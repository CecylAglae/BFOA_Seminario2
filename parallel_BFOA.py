from copy import deepcopy
from multiprocessing import Manager, Pool
import time
from BFOAMejorado import BFOA_mejorado
import numpy
import random
from fastaReader import fastaReader

if __name__ == "__main__":
    # Parámetros configurables
    numeroDeBacterias = 10       # Tamaño de la población
    numRandomBacteria = 2        # Bacterias que mantienen aleatoriedad
    iteraciones = 30             # Número de iteraciones
    tumbo = 2                    # Número de gaps a insertar en cada tumbo
    nado = 3                     # Pasos de nado después de tumbo
    k_orbc = 5                   # Frecuencia de aplicación del ORBC
    umbral_orbc = 0.7            # Umbral de sensibilidad al consenso
    prob_reemplazo_orbc = 0.3    # Probabilidad de aplicar ORBC a cada bacteria
    
    # Parámetros de interacción
    dAttr = 0.1                  # Coeficiente de atracción
    wAttr = 0.002                # Peso de atracción
    hRep = dAttr                 # Coeficiente de repulsión
    wRep = 0.001                 # Peso de repulsión

    # Lectura de secuencias
    fr = fastaReader()
    secuencias = fr.seqs
    names = fr.names
    
    # Convertir secuencias a listas de caracteres
    for i in range(len(secuencias)):
        secuencias[i] = list(secuencias[i])
        
    globalNFE = 0                # Contador de evaluaciones de función objetivo
    numSec = len(secuencias)     # Número de secuencias a alinear
    print(f"Número de secuencias: {numSec}")
    
    # Configuración de multiprocessing
    manager = Manager()
    poblacion = manager.list(range(numeroDeBacterias))
    names_manager = manager.list(names)
    NFE = manager.list([0]*numeroDeBacterias)

    # Función para inicializar población
    def poblacionInicial():
        for i in range(numeroDeBacterias):
            bacterium = []
            for j in range(numSec):
                bacterium.append(deepcopy(secuencias[j]))
            poblacion[i] = bacterium
   
    # Función auxiliar para mostrar población (debug)
    def printPoblacion():
        for i in range(numeroDeBacterias):
            print(f"Bacteria {i}:")
            for seq in poblacion[i]:
                print(''.join(seq))
            print("---")

    # Inicialización del operador bacterial mejorado
    operadorBacterial = BFOA_mejorado(numeroDeBacterias)
    veryBest = [None, None, None]  # [indice, fitness, secuencias]

    # Registro de tiempo de inicio
    start_time = time.time()
    
    print("Inicializando población...")
    poblacionInicial() 
    
    print("Iniciando proceso de optimización...")
    for it in range(iteraciones):
        print(f"\n--- Iteración {it+1}/{iteraciones} ---")
        
        # 1. Fase de tumbo (movimiento aleatorio)
        print("Aplicando tumbo...")
        operadorBacterial.tumbo(numSec, poblacion, tumbo)
        
        # 2. Cuadrar secuencias
        print("Cuadrando secuencias...")
        operadorBacterial.cuadra(numSec, poblacion)
        
        # 3. Aplicar operador de refinamiento basado en consenso
        print("Aplicando refinamiento por consenso...")
        operadorBacterial.aplicar_orbc(poblacion, numSec, k=k_orbc, 
                                      umbral=umbral_orbc, 
                                      prob_reemplazo=prob_reemplazo_orbc)
        
        # 4. Crear gran lista de pares para evaluación
        print("Creando lista de pares...")
        operadorBacterial.creaGranListaPares(poblacion)
        
        # 5. Evaluación paralela con BLOSUM
        print("Evaluando con BLOSUM...")
        operadorBacterial.evaluaBlosum()  # Paralelo
        
        # 6. Crear tablas de atracción y repulsión en paralelo
        print("Calculando interacciones bacterianas...")
        operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, hRep, wRep)
        
        # 7. Calcular interacción total y fitness
        print("Calculando fitness...")
        operadorBacterial.creaTablaInteraction()
        operadorBacterial.creaTablaFitness()
        
        # 8. Actualizar contador de evaluaciones
        globalNFE += operadorBacterial.getNFE()
        
        # 9. Identificar y almacenar mejor solución
        bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
        if (veryBest[0] is None) or (bestFitness > veryBest[1]):
            veryBest[0] = bestIdx
            veryBest[1] = bestFitness
            veryBest[2] = deepcopy(poblacion[bestIdx])
        
        # 10. Reemplazar peores soluciones
        operadorBacterial.replaceWorst(poblacion, veryBest[0])
        
        # 11. Reiniciar listas para siguiente iteración
        operadorBacterial.resetListas(numeroDeBacterias)
        
        # Mostrar progreso cada 5 iteraciones
        if (it+1) % 5 == 0:
            print(f"\nProgreso: {it+1}/{iteraciones} iteraciones")
            print(f"Mejor fitness actual: {veryBest[1]}")
            print(f"NFE acumuladas: {globalNFE}")
            print(f"Tiempo transcurrido: {time.time() - start_time:.2f} segundos")

    # Resultados finales
    print("\n--- Resultados Finales ---")
    print(f"Mejor solución encontrada (Fitness: {veryBest[1]}):")
    
    # Mostrar alineamiento final
    print("\nAlineamiento óptimo:")
    for i, seq in enumerate(veryBest[2]):
        print(f"{names[i]}: {''.join(seq)}")
    
    # Métricas de rendimiento
    print("\n--- Métricas de Rendimiento ---")
    print(f"Tiempo total de ejecución: {time.time() - start_time:.2f} segundos")
    print(f"Evaluaciones de función (NFE): {globalNFE}")
    print(f"Fitness alcanzado: {veryBest[1]}")
    
    # Guardar resultados en archivo
    with open("mejor_alineamiento.txt", "w") as f:
        f.write(f"Fitness: {veryBest[1]}\n\n")
        for i, seq in enumerate(veryBest[2]):
            f.write(f">{names[i]}\n{''.join(seq)}\n")
    
    print("\nResultados guardados en 'mejor_alineamiento.txt'")