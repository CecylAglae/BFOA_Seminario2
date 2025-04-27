from copy import deepcopy
import random
from bacteria import bacteria
from evaluadorBlosum import evaluadorBlosum

class BFOA_mejorado(bacteria):
    def __init__(self, numBacterias):
        super().__init__(numBacterias)
        self.mejor_consenso = None
        self.contador_iter = 0
    
    def aplicar_orbc(self, poblacion, numSec, k=5, umbral=0.7, prob_reemplazo=0.4):
        self.contador_iter += 1
        
        if self.contador_iter % k == 0 and len(poblacion) > 1 and numSec > 0:
            try:
                # Paso 1: Construir perfil de consenso
                consenso = self.construir_perfil_consenso(poblacion, numSec)
                if not consenso:
                    return
                
                # Paso 2: Alinear bacterias seleccionadas al consenso
                for i in range(len(poblacion)):
                    if random.random() < prob_reemplazo:
                        poblacion[i] = self.realinear_por_consenso(poblacion[i], consenso, umbral)
                
                # Actualizar mejor consenso histórico
                current_score = self.evaluar_consenso(consenso)
                if (self.mejor_consenso is None) or (current_score > self.evaluar_consenso(self.mejor_consenso)):
                    self.mejor_consenso = consenso
            except Exception as e:
                print(f"Error en ORBC: {str(e)}")
                return
    
    def construir_perfil_consenso(self, poblacion, num_secuencias):
        # Verificar población válida
        if not poblacion or num_secuencias <= 0:
            return []
        
        # Obtener longitud mínima segura
        longitudes = []
        for bacteria in poblacion:
            if bacteria and len(bacteria) >= num_secuencias:
                for sec in range(num_secuencias):
                    if sec < len(bacteria):
                        longitudes.append(len(bacteria[sec]))
        
        if not longitudes:
            return []
        
        longitud = min(longitudes)
        
        # Inicializar matriz de frecuencias
        perfil = {}
        for pos in range(longitud):
            perfil[pos] = {'A':0, 'T':0, 'C':0, 'G':0, '-':0}
            for bacteria in poblacion:
                if bacteria and len(bacteria) >= num_secuencias:
                    for sec in range(num_secuencias):
                        if sec < len(bacteria) and pos < len(bacteria[sec]):
                            nucleotido = bacteria[sec][pos]
                            perfil[pos][nucleotido] += 1
        
        # Determinar consenso por posición
        consenso = []
        for pos in sorted(perfil.keys()):
            max_nuc = max(perfil[pos].items(), key=lambda x: x[1])[0]
            consenso.append(max_nuc)
        
        return consenso
    
    def realinear_por_consenso(self, bacteria, consenso, umbral):
        if not bacteria or not consenso:
            return bacteria
        
        nueva_bacteria = deepcopy(bacteria)
        cons_len = len(consenso)
        
        for i in range(len(nueva_bacteria)):
            if i >= len(nueva_bacteria):
                continue
                
            secuencia = nueva_bacteria[i]
            nueva_secuencia = []
            
            for pos in range(cons_len):
                if pos < len(secuencia) and pos < cons_len:
                    if secuencia[pos] == consenso[pos]:
                        nueva_secuencia.append(secuencia[pos])
                    else:
                        if random.random() > umbral:
                            nueva_secuencia.append('-')
                        else:
                            nueva_secuencia.append(secuencia[pos])
                else:
                    nueva_secuencia.append('-')
            
            # Asegurar longitud correcta
            nueva_secuencia = nueva_secuencia[:cons_len]
            if len(nueva_secuencia) < cons_len:
                nueva_secuencia.extend(['-']*(cons_len - len(nueva_secuencia)))
            
            nueva_bacteria[i] = nueva_secuencia
        
        return nueva_bacteria
    
    def evaluar_consenso(self, consenso):
        if not consenso:
            return float('-inf')
            
        evaluador = evaluadorBlosum()
        score = 0
        valid_pairs = 0
        
        for i in range(len(consenso)-1):
            if i+1 < len(consenso):
                score += evaluador.getScore(consenso[i], consenso[i+1])
                valid_pairs += 1
        
        return score / valid_pairs if valid_pairs > 0 else 0