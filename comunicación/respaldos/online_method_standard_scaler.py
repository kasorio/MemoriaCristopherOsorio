# --------------------------------------------------------------------------
# 1. CONFIGURACIÓN INICIAL: SILENCIAR LOGS Y ADVERTENCIAS
# --------------------------------------------------------------------------
import os
import warnings
import sys
import traceback

# Silenciar logs de bajo nivel (C++) de TensorFlow
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
warnings.filterwarnings('ignore')

# --------------------------------------------------------------------------
# 2. IMPORTACIÓN DE LIBRERÍAS
# --------------------------------------------------------------------------
import json
import numpy as np
import time
from keras.models import load_model
import joblib  # NUEVO: Importar joblib para cargar el scaler

# --------------------------------------------------------------------------
# 3. CARGA DEL MODELO Y EL SCALER (ADAPTADO PARA UN SOLO HILO)
# --------------------------------------------------------------------------
# print("--- Iniciando el Servidor de Inferencia con Tuberías (Pipes) ---")
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # --- Carga del Modelo ---
    #MODEL_FILENAME = 'mlp_model3_500sim_antigua_f(x)_prob_250ep_1000bs_5co_relu.h5'
    MODEL_FILENAME = 'mlp_model3_symlog_24_11_25.h5'
    MODEL_PATH = os.path.join(script_dir, MODEL_FILENAME)
    # print(f"Cargando modelo desde: {MODEL_PATH}")
    model = load_model(MODEL_PATH)
    # print("Modelo cargado exitosamente.")

    # --- NUEVO: Carga del Scaler ---
    # ¡IMPORTANTE! Asegúrate de que este nombre de archivo coincida con tu archivo scaler.
    #SCALER_FILENAME = 'scaler_model3_500sim_antigua_f(x)_prob_250ep_1000bs_5co_relu.gz' 
    SCALER_FILENAME = 'scaler_model3_symlog_24_11_25.gz'
    SCALER_PATH = os.path.join(script_dir, SCALER_FILENAME)
    # print(f"Cargando escalador desde: {SCALER_PATH}")
    scaler = joblib.load(SCALER_PATH)
    # print("Escalador cargado exitosamente.")

except Exception as e:
    print(f'FATAL: Error al cargar los artefactos del modelo/scaler: {e}', file=sys.stderr)
    traceback.print_exc()
    sys.exit(1)

# --------------------------------------------------------------------------
#   FUNCION AUXILIAR PARA TRANSFORMAR Y TRATAR RANGOS DE DATOS DE MEJOR FORMA
# --------------------------------------------------------------------------
def symlog_transform(x):
    """Logaritmo simétrico: sign(x) * log(1 + |x|)"""
    # Nota: Usamos np.abs y np.log1p para seguridad numérica
    return np.sign(x) * np.log1p(np.abs(x))

# --------------------------------------------------------------------------
# 4. FUNCIÓN PARA PROCESAR PREDICCIONES (MODIFICADA CON SCALER)
# --------------------------------------------------------------------------
def process_prediction(input_json):
    try:
        data = json.loads(input_json)
        features_vector = data['features']
        
        if not isinstance(features_vector, list) or len(features_vector) != 7:
            raise ValueError(f"Se esperaba una lista de 7 características")
        
        # 1. Convertir a array de NumPy (Shape: 1x7)
        features = np.array(features_vector, dtype=np.float64).reshape(1, -1)

        # 2. Validar integridad básica
        if np.any(np.isnan(features)) or np.any(np.isinf(features)):
            # Opcional: Podrías intentar reparar infinitos aquí si C++ los manda,
            # pero idealmente C++ manda números grandes finitos.
            raise ValueError("El vector de características contiene NaN o Inf")

        # 3. --- APLICAR SYMLOG TRANSFORM ---
        # Aplicamos SOLO a las columnas que transformamos en el entrenamiento.
        # Índices: 0:vars, 1:restr, 2:lb, 3:ub, 4:space, 5:big_d, 6:low_d
        indices_to_transform = [2, 3, 4, 5, 6]
        
        for i in indices_to_transform:
            features[0, i] = symlog_transform(features[0, i])

        # 4. --- APLICAR SCALER ---
        scaled_features = scaler.transform(features)
        
        # 5. Predicción
        probabilities = model.predict(scaled_features, verbose=0)[0]

        decision_index = int(np.argmax(probabilities))
        
        print(f"bis: {decision_index}, id: {data.get('id', -1)}")

        response = {
            'decision': decision_index,
            'id': data.get('id', -1),
            'probabilities': probabilities.tolist()
        }
        
        return json.dumps(response)
                
    except (json.JSONDecodeError, KeyError, ValueError) as e:
        error_id = -1
        try:
            data = json.loads(input_json)
            error_id = data.get('id', -1)
        except:
            pass
            
        print(f"Error de procesamiento para ID {error_id}: {e}")
        response = {
            'error': f'Error de procesamiento de datos: {str(e)}',
            'id': error_id,
            'timestamp': time.time()
        }
        return json.dumps(response)
        
    except Exception as e:
        error_id = -1
        try:
            data = json.loads(input_json)
            error_id = data.get('id', -1)
        except:
            pass
        
        print(f"Error inesperado para ID {error_id}: {str(e)}")
        traceback.print_exc()
        response = {
            'error': f'Error inesperado: {str(e)}',
            'id': error_id,
            'timestamp': time.time()
        }
        return json.dumps(response)


# --------------------------------------------------------------------------
# 5. SERVIDOR PRINCIPAL CON TUBERÍAS (SIN CAMBIOS)
# --------------------------------------------------------------------------
def start_pipe_server():
    """
    Inicia un servidor que escucha peticiones en una tubería con nombre (FIFO)
    y envía respuestas a través de otra.
    """
    REQUEST_PIPE_PATH = "/tmp/ibex_request_pipe"
    RESPONSE_PIPE_PATH = "/tmp/ibex_response_pipe"
    
    # Asegurarse de que las tuberías existen
    if not os.path.exists(REQUEST_PIPE_PATH):
        print(f"Error: La tubería de peticiones no existe en '{REQUEST_PIPE_PATH}'", file=sys.stderr)
        print("Por favor, créala con: mkfifo /tmp/ibex_request_pipe", file=sys.stderr)
        sys.exit(1)
        
    if not os.path.exists(RESPONSE_PIPE_PATH):
        print(f"Error: La tubería de respuestas no existe en '{RESPONSE_PIPE_PATH}'", file=sys.stderr)
        print("Por favor, créala con: mkfifo /tmp/ibex_response_pipe", file=sys.stderr)
        sys.exit(1)

    # print(f"Servidor escuchando en la tubería: {REQUEST_PIPE_PATH}")
    # print("Presione Ctrl+C para detener el servidor.")

    try:
        while True:
            with open(REQUEST_PIPE_PATH, 'r') as request_pipe:
                request_line = request_pipe.readline()

            if request_line:
                response_json = process_prediction(request_line.strip())
                with open(RESPONSE_PIPE_PATH, 'w') as response_pipe:
                    response_pipe.write(response_json + '\n')
    
    except KeyboardInterrupt:
        print("\nCerrando servidor por solicitud del usuario...")
    except Exception as e:
        print(f"Error crítico en el servidor de tuberías: {e}", file=sys.stderr)
        traceback.print_exc()
    finally:
        print("Servidor detenido.")

# --------------------------------------------------------------------------
# 6. FUNCIÓN PRINCIPAL (SIN CAMBIOS)
# --------------------------------------------------------------------------
if __name__ == '__main__':
    start_pipe_server()