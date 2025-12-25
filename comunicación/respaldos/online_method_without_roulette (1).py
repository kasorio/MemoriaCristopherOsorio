# --------------------------------------------------------------------------
# 1. CONFIGURACIÓN INICIAL: SILENCIAR LOGS Y ADVERTENCIAS
# --------------------------------------------------------------------------
import os
import warnings
import sys

# Silenciar logs de bajo nivel (C++) de TensorFlow
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
# Opcional: Forzar el uso de la CPU
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
# Suprimir advertencias de Python
warnings.filterwarnings('ignore')

# --------------------------------------------------------------------------
# 2. IMPORTACIÓN DE LIBRERÍAS
# --------------------------------------------------------------------------
import json
import numpy as np
import joblib
import socket
import threading
from keras.models import load_model
from concurrent.futures import ThreadPoolExecutor

MAX_WORKERS = min(32, (os.cpu_count() or 1) + 1)

# --------------------------------------------------------------------------
# 3. CARGA DEL MODELO UNA SOLA VEZ
# --------------------------------------------------------------------------
try:
    # Ruta del directorio del script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    MODEL_PATH = os.path.join(script_dir, 'mlp_model3_500sim_without_scaler.h5')
#SCALER_PATH = os.path.join(script_dir, 'scaler_model3_500sim.gz')
    
    #print(f"Cargando modelo desde: {MODEL_PATH}")
    #print(f"Cargando scaler desde: {SCALER_PATH}")
    
    model = load_model(MODEL_PATH)
    #scaler = joblib.load(SCALER_PATH)

    # print(f"media del scaler: {scaler.mean_}")
    # print(f"std del scaler: {scaler.scale_}")
    
    #print("Modelo y scaler cargados exitosamente")
    
except Exception as e:
    print(f'Error al cargar el modelo: {e}', file=sys.stderr)
    sys.exit(1)

# --------------------------------------------------------------------------
# 4. FUNCIÓN PARA PROCESAR PREDICCIONES
# --------------------------------------------------------------------------
def process_prediction(input_json):
    """
    Procesa una predicción basada en el JSON de entrada
    """
    try:
        data = json.loads(input_json)
        features_vector = data['features']

        #print(f"Datos recibidos: {features_vector}")
        
        features = np.array(features_vector).reshape(1, -1)

        #print(f"Características: {features}")

        # Verificar que el vector tiene 7 características
        if len(features_vector) != 7:
            raise ValueError(f"Se esperaba un vector de 7 características, pero llegaron {len(features_vector)}")
        
        # Realizar la predicción
        # La salida es un arreglo de probabilidades [LS, LF, RR, SM, SS, SSR]
        probabilities = model.predict(features, verbose=0)[0]

        # Elegir la clase con mayor probabilidad
        decision_index = int(np.argmax(probabilities))

        print(f"Probabilidades: {probabilities}")
        print(f"response: {decision_index}")

        # Preparar la respuesta
        response = {
            'decision': decision_index,
            'id': data.get('id', -1),
            'probabilities': probabilities.tolist()  # Opcional: incluir probabilidades
        }
        #print(f"Respuesta: {response}")
        return json.dumps(response)    
        
    except (json.JSONDecodeError, KeyError, ValueError) as e:
        error_id = -1
        try:
            data = json.loads(input_json)
            error_id = data.get('id', -1)
        except:
            pass
            
        response = {
            'error': f'Error de procesamiento de datos: {e}',
            'id': error_id
        }
        return json.dumps(response)
        
    except Exception as e:
        error_id = -1
        try:
            data = json.loads(input_json)
            error_id = data.get('id', -1)
        except:
            pass
            
        response = {
            'error': f'Error inesperado: {str(e)}',
            'id': error_id
        }
        return json.dumps(response)

# --------------------------------------------------------------------------
# 5. MANEJADOR DE CLIENTE
# --------------------------------------------------------------------------
def handle_client(client_socket, address):
    """
    Maneja la conexión de un cliente
    """
    #print(f"Conexión establecida desde {address}")
    
    try:
        while True:
            # Recibir datos del cliente
            data = client_socket.recv(1024)
            if not data:
                break
                
            # Decodificar el JSON
            input_json = data.decode('utf-8')
            #print(f"Datos recibidos de {address}: {input_json}")
            
            # Procesar la predicción
            response = process_prediction(input_json)
            #print(f"Respuesta enviada a {address}: {response}")
            
            # Enviar respuesta
            client_socket.send(response.encode('utf-8'))
            
    except Exception as e:
        print(f"Error manejando cliente {address}: {e}")
        error_response = json.dumps({
            'error': f'Error del servidor: {str(e)}',
            'id': -1
        })
        try:
            client_socket.send(error_response.encode('utf-8'))
        except:
            pass
    finally:
        client_socket.close()
        #print(f"Conexión cerrada con {address}")

# --------------------------------------------------------------------------
# 6. SERVIDOR PRINCIPAL
# --------------------------------------------------------------------------
def start_server(host='localhost', port=8888):
    """
    Inicia el servidor socket con ThreadPoolExecutor
    """
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    
    # NUEVO: Usar ThreadPoolExecutor para mejor gestión
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        try:
            server_socket.bind((host, port))
            server_socket.listen(5)
            #print(f"Servidor iniciado en {host}:{port}")
            #print(f"Usando ThreadPoolExecutor con {MAX_WORKERS} workers máximo")
            #print("Esperando conexiones...")
            
            while True:
                client_socket, address = server_socket.accept()
                # NUEVO: Usar executor.submit en lugar de threading.Thread
                executor.submit(handle_client, client_socket, address)
                
        except KeyboardInterrupt:
            print("\nCerrando servidor...")
        except Exception as e:
            print(f"Error en el servidor: {e}")
        finally:
            server_socket.close()

# --------------------------------------------------------------------------
# 7. FUNCIÓN PRINCIPAL
# --------------------------------------------------------------------------
def main():
    """
    Función principal del servidor
    """
    # Configuración del servidor
    HOST = 'localhost'
    PORT = 8888
    
    # Verificar argumentos de línea de comandos
    if len(sys.argv) > 1:
        try:
            PORT = int(sys.argv[1])
        except ValueError:
            print("Puerto inválido. Usando puerto por defecto 8888")
    
    if len(sys.argv) > 2:
        HOST = sys.argv[2]
    
    #print(f"Iniciando servidor de modelo ML en {HOST}:{PORT}")
    start_server(HOST, PORT)

if __name__ == '__main__':
    main()