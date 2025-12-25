//                                  I B E X
// File        : ibex_Optimizer.cpp
// Author      : Gilles Chabert, Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : Feb 13, 2025
//============================================================================

#include "ibex_Optimizer.h"
#include "ibex_Timer.h"
#include "ibex_Function.h"
#include "ibex_NoBisectableVariableException.h"
#include "ibex_BxpOptimData.h"
#include "ibex_CovOptimData.h"
#include "ibex_CellBeamSearch.h"
#include "ibex_LoupFinderDefault.h"

#include "ibex_SmearFunction.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_OptimLargestFirst.h"
#include "ibex_System.h"

#include <float.h>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <iomanip>
#include <chrono>

#include <limits> // Necesario para std::numeric_limits
#include <cmath>  // Necesario para std::abs

using namespace std;

// ============================================================================
//         COMUNICACIÓN CON MODELO MEDIANTE SOCKET
// ============================================================================
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <string>
#include <iostream>

// Variable global para mantener la conexión
static int model_socket = -1;
static bool socket_connected = false;

bool connect_to_model_server(const std::string& host, int port) {
    if (socket_connected) return true;
    
    model_socket = socket(AF_INET, SOCK_STREAM, 0);
    if (model_socket == -1) {
        std::cerr << "Error creating socket" << std::endl;
        return false;
    }
    
    struct sockaddr_in serv_addr;
    memset(&serv_addr, 0, sizeof(serv_addr));  // Initialize structure
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(port);
    
    if (inet_pton(AF_INET, host.c_str(), &serv_addr.sin_addr) <= 0) {
        std::cerr << "Invalid address: " << host << std::endl;
        close(model_socket);
        model_socket = -1;
        return false;
    }
    
    if (connect(model_socket, (struct sockaddr*)&serv_addr, sizeof(serv_addr)) < 0) {
        std::cerr << "Connection failed to " << host << ":" << port << std::endl;
        close(model_socket);
        model_socket = -1;
        return false;
    }
    
    socket_connected = true;
    //std::cout << "Connected to model server at " << host << ":" << port << std::endl;
    return true;
}

// Overloaded function with default parameters
bool connect_to_model_server() {
    return connect_to_model_server("localhost", 8888);
}

// Búfer estático para almacenar datos de respuesta incompletos entre llamadas
static std::string response_buffer;

std::string call_python_model(const std::string& input_json) {
    // Intentar conectar si no estamos conectados
    if (!socket_connected && !connect_to_model_server()) {
        return "{\"error\": \"No connection to model server\"}";
    }
    
    // --- 1. MODIFICACIÓN AL ENVIAR ---
    // Añadimos el delimitador de nueva línea '\n' al final del JSON.
    std::string message_to_send = input_json + "\n";
    
    ssize_t bytes_sent = send(model_socket, message_to_send.c_str(), message_to_send.length(), 0);
    if (bytes_sent == -1) {
        std::cerr << "Failed to send data" << std::endl;
        socket_connected = false;
        close(model_socket);
        model_socket = -1;
        return "{\"error\": \"Failed to send data\"}";
    }
    
    // --- 2. MODIFICACIÓN AL RECIBIR ---
    // Usamos un bucle para leer del socket hasta encontrar un mensaje completo ('\n').
    while (true) {
        // Primero, revisamos si ya tenemos un mensaje completo en nuestro búfer
        size_t newline_pos = response_buffer.find('\n');
        if (newline_pos != std::string::npos) {
            // ¡Sí! Tenemos un mensaje completo.
            std::string full_response = response_buffer.substr(0, newline_pos);
            // Eliminamos el mensaje y el '\n' del búfer.
            response_buffer.erase(0, newline_pos + 1);
            return full_response;
        }
        
        // Si no, necesitamos leer más datos del socket.
        char temp_buffer[1024] = {0};
        int bytes_received = recv(model_socket, temp_buffer, sizeof(temp_buffer) - 1, 0);
        
        if (bytes_received > 0) {
            // Añadimos los nuevos datos a nuestro búfer de respuesta.
            response_buffer.append(temp_buffer, bytes_received);
        } else {
            // Hubo un error o el servidor cerró la conexión.
            std::cerr << "Failed to receive response or connection closed" << std::endl;
            socket_connected = false;
            close(model_socket);
            model_socket = -1;
            response_buffer.clear(); // Limpiamos el búfer en caso de error
            return "{\"error\": \"Failed to receive response\"}";
        }
    }
}

void close_model_connection() {
    if (socket_connected && model_socket != -1) {
        close(model_socket);
        socket_connected = false;
        model_socket = -1;
        //std::cout << "Model server connection closed" << std::endl;
    }
}
// ============================================================================
//                      FIN DE LA COMUNICACIÓN CON EL MODELO
// ============================================================================

// ============================================================================
//  AQUÍ SE REALIZA EL PREPROCESAMIENTO DE LOS DATOS ANTES DE LLAMAR AL MODELO
// ============================================================================

double const UMBRAL_CERO = 1e-7;

enum class FeatureType {
    LbFObj, UbFObj, SearchSpace, BiggerDiam, LowerDiam, Loup,
    GapRelLoup, BufferSize, Depth, Variables, Restricciones,
    GenericFeature
};

double preprocess_feature(double val, FeatureType type) {
    // 1. Manejar NaN. Reemplazar con 0 puede ocultar problemas. 
    //    Considera si otro valor (como -1) o un manejo de errores es mejor.
    if (std::isnan(val)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // 2. Tratar valores muy cercanos a cero
    if (std::abs(val) < UMBRAL_CERO) {
        val = 0.0;
    }

    switch (type) {
		case FeatureType::LbFObj:
			// Límite inferior y superior (-1e7, 1e7)
			return std::max(-10000000.0, std::min(val, 10000000.0));

        case FeatureType::UbFObj:
            return std::min(val, 1000000.0);

        // case FeatureType::SearchSpace:
        //     // Límite superior de 500000.0
        //     return std::min(val, 500000.0);
            
        case FeatureType::BiggerDiam:
            // Límite superior de 500000.0
            return std::min(val, 1000000.0);

		// case FeatureType::LowerDiam:
		// 	// Límite inferior de -500000.0
		// 	return std::min(val, 10000000.0);

        // case FeatureType::Loup:
		// 	break;
    
		// case FeatureType::BufferSize:
        //     break;
    
        default:
            // Para 'gap_rel_loup', 'depth', 'variables', etc.
            // Si el valor es infinito (y no fue "clippeado" arriba), es un problema.
            // Lo reemplazamos por un valor muy grande pero finito.
            if (std::isinf(val)) {
                return (val > 0) ? std::numeric_limits<double>::max() : -std::numeric_limits<double>::max();
            }
            break;
    }

    // Si no es un caso especial con clipping, devolver el valor después del chequeo de cero.
    return val;
}

FeatureType string_to_feature_type(const std::string& feature_name) {
    static const std::unordered_map<std::string, FeatureType> map = {
        {"lb_f_obj",     FeatureType::LbFObj},
        {"ub_f_obj",     FeatureType::UbFObj},
        // {"search_space", FeatureType::SearchSpace},
        {"bigger_diam",  FeatureType::BiggerDiam},
        // {"lower_diam",   FeatureType::LowerDiam}
        // Añade aquí los mismos nombres de clave que usas en Python
    };

    auto it = map.find(feature_name);
    if (it != map.end()) {
        return it->second;
    }
    return FeatureType::GenericFeature;
}

double preprocess_feature(double raw_value, const std::string& feature_name) {
    // Esta función ahora sabe que existe otra 'preprocess_feature' que acepta un 'FeatureType'.
    // El compilador elegirá la correcta y no intentará una conversión inválida.
    return preprocess_feature(raw_value, string_to_feature_type(feature_name));
}

// ================================================================================
//  FIN AQUÍ SE REALIZA EL PREPROCESAMIENTO DE LOS DATOS ANTES DE LLAMAR AL MODELO
// ================================================================================


namespace ibex {


/*
 * TODO: redundant with ExtendedSystem.
 */
void Optimizer::write_ext_box(const IntervalVector& box, IntervalVector& ext_box) {
	int i2=0;
	for (int i=0; i<n; i++,i2++) {
		if (i2==goal_var) i2++; // skip goal variable
		ext_box[i2]=box[i];
	}
}

void Optimizer::read_ext_box(const IntervalVector& ext_box, IntervalVector& box) {
	int i2=0;
	for (int i=0; i<n; i++,i2++) {
		if (i2==goal_var) i2++; // skip goal variable
		box[i]=ext_box[i2];
	}
}

Optimizer::Optimizer(int n, Ctc& ctc, Bsc& bsc, LoupFinder& finder,
		CellBufferOptim& buffer,
		int goal_var, double eps_x, double rel_eps_f, double abs_eps_f,
		bool enable_statistics) :
                						n(n), goal_var(goal_var),
										ctc(ctc), bsc(bsc), loup_finder(finder), buffer(buffer),
										eps_x(n, eps_x), rel_eps_f(rel_eps_f), abs_eps_f(abs_eps_f),
										trace(0), timeout(-1), extended_COV(true), anticipated_upper_bounding(true),
										status(SUCCESS),
										uplo(NEG_INFINITY), uplo_of_epsboxes(POS_INFINITY), loup(POS_INFINITY),
										loup_point(IntervalVector::empty(n)), initial_loup(POS_INFINITY), loup_changed(false),
										time(0), nb_cells(0), cov(NULL), python_time(0), wall_time(0){

	if (trace) std::cout.precision(12);
	
	if (enable_statistics) {
		statistics = new Statistics();
		// TODO: enable statistics for missing operators (cell buffer)
		bsc.enable_statistics(*statistics, "Bsc"); 
		ctc.enable_statistics(*statistics, "Ctc"); 
		loup_finder.enable_statistics(*statistics, "LoupFinder"); 
	} else
		statistics = NULL;
}

Optimizer::Optimizer(OptimizerConfig& config) :
	Optimizer(
		config.nb_var(), 
		config.get_ctc(), 
		config.get_bsc(), 
		config.get_loup_finder(),
		config.get_cell_buffer(),
		config.goal_var(),
		OptimizerConfig::default_eps_x, // tmp, see below
		config.get_rel_eps_f(),
		config.get_abs_eps_f(),
		config.with_statistics()) {

	(Vector&) eps_x				= config.get_eps_x();
	trace						= config.get_trace();
	timeout						= config.get_timeout();
	extended_COV				= config.with_extended_cov();
	anticipated_upper_bounding	= config.with_anticipated_upper_bounding();
}

Optimizer::~Optimizer() {
	if (cov) delete cov;
	if (statistics) delete statistics;
}

// compute the value ymax (decreasing the loup with the precision)
// the heap and the current box are contracted with y <= ymax
double Optimizer::compute_ymax() {
	if (anticipated_upper_bounding) {
		//double ymax = loup - rel_eps_f*fabs(loup); ---> wrong :the relative precision must be correct for ymax (not loup)
		double ymax = loup>0 ?
				1/(1+rel_eps_f)*loup
		:
				1/(1-rel_eps_f)*loup;

		if (loup - abs_eps_f < ymax)
			ymax = loup - abs_eps_f;
		//return ymax;
		return next_float(ymax);
	} else
		return loup;
}

bool Optimizer::update_loup(const IntervalVector& box, BoxProperties& prop) {

	try {

		pair<IntervalVector,double> p=loup_finder.find(box,loup_point,loup,prop);
		loup_point = p.first;
		loup = p.second;

		if (trace) {
			std::cout << "                    ";
			std::cout << "\033[32m loup= " << loup << "\033[0m" << endl;
//			cout << " loup point=";
//			if (loup_finder.rigorous())
//				cout << loup_point << endl;
//			else
//				cout << loup_point.lb() << endl;
		}
		return true;

	} catch(LoupFinder::NotFound&) {
		return false;
	}
}

//bool Optimizer::update_entailed_ctr(const IntervalVector& box) {
//	for (int j=0; j<m; j++) {
//		if (entailed->normalized(j)) {
//			continue;
//		}
//		Interval y=sys.ctrs[j].f.eval(box);
//		if (y.lb()>0) return false;
//		else if (y.ub()<=0) {
//			entailed->set_normalized_entailed(j);
//		}
//	}
//	return true;
//}

void Optimizer::update_uplo() {
	double new_uplo=POS_INFINITY;

	if (! buffer.empty()) {
		new_uplo= buffer.minimum();
		if (new_uplo > loup && uplo_of_epsboxes > loup) {
			std::cout << " loup = " << loup << " new_uplo=" << new_uplo <<  " uplo_of_epsboxes=" << uplo_of_epsboxes << endl;
			ibex_error("optimizer: new_uplo>loup (please report bug)");
		}
		if (new_uplo < uplo) {
			std::cout << "uplo= " << uplo << " new_uplo=" << new_uplo << endl;
			ibex_error("optimizer: new_uplo<uplo (please report bug)");
		}

		// uplo <- max(uplo, min(new_uplo, uplo_of_epsboxes))
		if (new_uplo < uplo_of_epsboxes) {
			if (new_uplo > uplo) {
				uplo = new_uplo;

				if (trace)
					std::cout << "\033[33m uplo= " << uplo << "\033[0m" << endl;
			}
		}
		else uplo = uplo_of_epsboxes;
	}
	else if (buffer.empty() && loup != POS_INFINITY) {
		// empty buffer : new uplo is set to ymax (loup - precision) if a loup has been found
		new_uplo=compute_ymax(); // not new_uplo=loup, because constraint y <= ymax was enforced
		//    cout << " new uplo buffer empty " << new_uplo << " uplo " << uplo << endl;

		double m = (new_uplo < uplo_of_epsboxes) ? new_uplo :  uplo_of_epsboxes;
		if (uplo < m) uplo = m; // warning: hides the field "m" of the class
		// note: we always have uplo <= uplo_of_epsboxes but we may have uplo > new_uplo, because
		// ymax is strictly lower than the loup.
	}

}

void Optimizer::update_uplo_of_epsboxes(double ymin) {

	// the current box cannot be bisected.  ymin is a lower bound of the objective on this box
	// uplo of epsboxes can only go down, but not under uplo : it is an upperbound for uplo,
	// that indicates a lowerbound for the objective in all the small boxes
	// found by the precision criterion
	assert (uplo_of_epsboxes >= uplo);
	assert(ymin >= uplo);
	if (uplo_of_epsboxes > ymin) {
		uplo_of_epsboxes = ymin;
		if (trace) {
			std::cout << " unprocessable tiny box: now uplo<=" << setprecision(12) <<  uplo_of_epsboxes << " uplo=" << uplo << endl;
		}
	}
}

void Optimizer::handle_cell(Cell& c) {

	contract_and_bound(c);

	if (c.box.is_empty()) {
		delete &c;
	} else {
		buffer.push(&c);
	}
}

void Optimizer::contract_and_bound(Cell& c) {

	/*======================== contract y with y<=loup ========================*/
	Interval& y=c.box[goal_var];

	double ymax;
	if (loup==POS_INFINITY) ymax = POS_INFINITY;
	// ymax is slightly increased to favour subboxes of the loup
	// TODO: useful with double heap??
	else ymax = compute_ymax()+1.e-15;

	y &= Interval(NEG_INFINITY,ymax);

	if (y.is_empty()) {
		c.box.set_empty();
		return;
	} else {
		c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
	}

	/*================ contract x with f(x)=y and g(x)<=0 ================*/
	//cout << " [contract]  x before=" << c.box << endl;
	//cout << " [contract]  y before=" << y << endl;

	ContractContext context(c.prop);
	if (c.bisected_var!=-1) {
		context.impact.clear();
		context.impact.add(c.bisected_var);
		context.impact.add(goal_var);
	}

	ctc.contract(c.box, context);
	//cout << c.prop << endl;
	if (c.box.is_empty()) return;

	//cout << " [contract]  x after=" << c.box << endl;
	//cout << " [contract]  y after=" << y << endl;
	/*====================================================================*/

	/*========================= update loup =============================*/

	IntervalVector tmp_box(n);
	read_ext_box(c.box,tmp_box);

	c.prop.update(BoxEvent(c.box,BoxEvent::CHANGE));

	bool loup_ch=update_loup(tmp_box, c.prop);

	// update of the upper bound of y in case of a new loup found
	if (loup_ch) {
		y &= Interval(NEG_INFINITY,compute_ymax());
		c.prop.update(BoxEvent(c.box,BoxEvent::CONTRACT,BitSet::singleton(n+1,goal_var)));
	}

	//TODO: should we propagate constraints again?

	loup_changed |= loup_ch;

	if (y.is_empty()) { // fix issue #44
		c.box.set_empty();
		return;
	}

	/*====================================================================*/
	// Note: there are three different cases of "epsilon" box,
	// - NoBisectableVariableException raised by the bisector (---> see optimize(...)) which
	//   is independent from the optimizer
	// - the width of the box is less than the precision given to the optimizer ("eps_x" for
	//   the original variables and "abs_eps_f" for the goal variable)
	// - the extended box has no bisectable domains (if eps_x=0 or <1 ulp)
	if (((tmp_box.diam()-eps_x).max()<=0 && y.diam() <=abs_eps_f) || !c.box.is_bisectable()) {
		update_uplo_of_epsboxes(y.lb());
		c.box.set_empty();
		return;
	}

	// ** important: ** must be done after upper-bounding
	//kkt.contract(tmp_box);

	if (tmp_box.is_empty()) {
		c.box.set_empty();
	} else {
		// the current extended box in the cell is updated
		write_ext_box(tmp_box,c.box);
	}
}

Optimizer::Status Optimizer::optimize(const IntervalVector& init_box, double obj_init_bound) {
	start(init_box, obj_init_bound);
	return optimize();
}


Optimizer::Status Optimizer::optimize(const CovOptimData& data, double obj_init_bound) {
	start(data, obj_init_bound);
	return optimize();
}

Optimizer::Status Optimizer::optimize(const char* cov_file, double obj_init_bound) {
	CovOptimData data(cov_file);
	start(data, obj_init_bound);
	return optimize();
}

void Optimizer::start(const IntervalVector& init_box, double obj_init_bound) {

	loup=obj_init_bound;

	// Just to initialize the "loup" for the buffer
	// TODO: replace with a set_loup function
	buffer.contract(loup);

	uplo=NEG_INFINITY;
	uplo_of_epsboxes=POS_INFINITY;

	nb_cells=0;

	buffer.flush();

	Cell* root=new Cell(IntervalVector(n+1));

	write_ext_box(init_box, root->box);

	// add data required by the bisector
	bsc.add_property(init_box, root->prop);

	// add data required by the contractor
	ctc.add_property(init_box, root->prop);

	// add data required by the buffer
	buffer.add_property(init_box, root->prop);

	// add data required by the loup finder
	loup_finder.add_property(init_box, root->prop);

	//cout << "**** Properties ****\n" << root->prop << endl;

	loup_changed=false;
	initial_loup=obj_init_bound;

	loup_point = init_box; //.set_empty();
	time=0;
	python_time = 0.0;
	wall_time = 0.0;
	start_wall_time = std::chrono::high_resolution_clock::now();

	// conectamos al servidor del modelo
	connect_to_model_server("127.0.0.1", 8888);

	if (cov) delete cov;
	cov = new CovOptimData(extended_COV? n+1 : n, extended_COV);
	cov->data->_optim_time = 0;
	cov->data->_optim_nb_cells = 0;

	handle_cell(*root);
}

void Optimizer::start(const CovOptimData& data, double obj_init_bound) {

	loup=obj_init_bound;

	// Just to initialize the "loup" for the buffer
	// TODO: replace with a set_loup function
	buffer.contract(loup);

	uplo=data.uplo();
	loup=data.loup();
	loup_point=data.loup_point();
	uplo_of_epsboxes=POS_INFINITY;

	nb_cells=0;

	buffer.flush();

	for (size_t i=loup_point.is_empty()? 0 : 1; i<data.size(); i++) {

		IntervalVector box(n+1);

		if (data.is_extended_space())
			box = data[i];
		else {
			write_ext_box(data[i], box);
			box[goal_var] = Interval(uplo,loup);
			ctc.contract(box);
			if (box.is_empty()) continue;
		}

		Cell* cell=new Cell(box);

		// add data required by the cell buffer
		buffer.add_property(box, cell->prop);

		// add data required by the bisector
		bsc.add_property(box, cell->prop);

		// add data required by the contractor
		ctc.add_property(box, cell->prop);

		// add data required by the loup finder
		loup_finder.add_property(box, cell->prop);

		buffer.push(cell);
	}

	loup_changed=false;
	initial_loup=obj_init_bound;

	time=0;
	python_time = 0.0;
	wall_time = 0.0;
	start_wall_time = std::chrono::high_resolution_clock::now();

	// conectamos al servidor del modelo
	connect_to_model_server("127.0.0.1", 8888);

	if (cov) delete cov;
	cov = new CovOptimData(extended_COV? n+1 : n, extended_COV);
	cov->data->_optim_time = data.time();
	cov->data->_optim_nb_cells = data.nb_cells();
}


/**
 * CALCULO DE LA FEATURE GAP REL LOUP
 * EL FEATURE ES MÁS ROBUSTO QUE EL ANTERIOR, POR LO QUE DEBIESE FUNCIONAR MEJOR
 * \brief Calculate a more robust relative gap feature based on the loup.
 */

double calculate_improved_gap_rel_loup(double current_loup, double lb_f_obj, double ub_f_obj) {
    double epsilon = 1e-12; // Valor más pequeño para mayor precisión
    double gap_rel_loup;
    
    if (!std::isfinite(current_loup) || current_loup == POS_INFINITY) {
        // Caso 1: No hay loup válido
        // El gap representa la incertidumbre relativa de la caja actual
        double numerator = ub_f_obj - lb_f_obj;
        double denominator = std::max(fabs(ub_f_obj), fabs(lb_f_obj));
		std::cout << " gap_rel_loup de ub - lb " << gap_rel_loup << std::endl;
        
        if (denominator < epsilon) {
            gap_rel_loup = 1.0; // Máxima incertidumbre cuando ambos límites son ~0
        } else {
			std::cout << " entro al else (denominator > epsilon) " << std::endl;
            gap_rel_loup = numerator / (denominator + epsilon);
        }
        
        // Normalizar a [0,1] - mayor valor = menos prometedora
        //gap_rel_loup = std::min(1.0, std::max(0.0, gap_rel_loup));
        
    } else {
        // Caso 2: Hay loup válido
        // Calculamos qué tan cerca está el lower bound de la mejor solución conocida
        double numerator = current_loup - lb_f_obj;
        
        // Si lb_f_obj > current_loup, la caja es muy prometedora (gap negativo)
        if (numerator < 0) {
            gap_rel_loup = 0.0; // Muy prometedora
        } else {
            // Normalizamos por la magnitud del loup
            double denominator = fabs(current_loup) + epsilon;
            gap_rel_loup = numerator / denominator;
			std::cout << " gap_rel_loup de loup - lb " << gap_rel_loup << std::endl;
            
            // Limitamos a [0,1] para evitar valores extremos
            //gap_rel_loup = std::min(1.0, gap_rel_loup);
        }
    }
    
    return gap_rel_loup;
}

Optimizer::Status Optimizer::optimize() {
	Timer timer;
	auto now = std::chrono::high_resolution_clock::now();
	timer.start();
	
	Timer python_timer;
	
	long call_id = 0;
	
	pid_t pid1 = getpid();

	update_uplo();

	try {
		/***********************************************
		 * DECLARACIÓN DE VARIABLES, BUFFERS Y FINDERS *
		 ***********************************************/
		CellBeamSearch * thebuffer = dynamic_cast<CellBeamSearch*>(&buffer);
		LoupFinderDefault * lfd = dynamic_cast<LoupFinderDefault*>(&loup_finder);

		// cola/vector que va a guardar las celdas
		queue<Cell*> aux;

		// variable aux para guardar el valor de la caja inicial
		IntervalVector aux_box = thebuffer->top()->box;
		double prec = 1e-7; //precisión para el solver (compara ub-lb y prec)
		/***************************************************
		 * FIN DECLARACIÓN DE VARIABLES, BUFFERS Y FINDERS *
		 ***************************************************/

		/*****************************
		 * DECLARACIÓN DE BISECTORES *
		 *****************************/
		OptimLargestFirst bisector_olf(goal_var, true, prec, 0.5);
		RoundRobin bisector_rr(prec, 0.5);

		/**********************************************************
		 * OBTENCIÓN DEL SISTEMA PARA BISECTORES QUE LO NECESITAN *
		 **********************************************************/
		System system = lfd->finder_x_taylor.sys;
		SmearMax bisector_sm(system, prec);
		SmearSum bisector_ss(system, prec);
		SmearSumRelative bisector_ssr(system, prec);
		/*********************************
	 	 * FIN DECLARACIÓN DE BISECTORES *
	 	 *********************************/
		
		/***************************
		 * BORDES Y REGIÓN INTERNA *
		 ***************************/ 
		double aux_uplo = uplo;
		double aux_loup = loup;
		IntervalVector inner = aux_box;
		/*******************************
		 * FIN BORDES Y REGIÓN INTERNA *
		 *******************************/

		//first node to process is pushed into the queue
		aux.push(thebuffer->top()); 

		
		// This variable is used to check the number of feasible divings realized.
		// It is used to call the MLP model
		int cont = 0;
		double epsilon = 1e-10;

		/**************************************************************
		 * ESTO ES LO QUE NO SE COMENTA, ES EL PROCEDIMIENTO COMO TAL *
		 **************************************************************/
		Cell* current_cell = aux.front();
		double current_loup = loup;
		double numerator;
		double denominator;
		
		//calculamos las features
		//int restricciones = lfd->finder_x_taylor.sys.nb_ctr;
		int variables = n; // cantidad de variables que tiene el problema de optimización
		double lb_f_obj = (thebuffer->top()->box)[goal_var].lb(); // límite inferior de la función objetivo evaluada en la caja actual
		double ub_f_obj = (thebuffer->top()->box)[goal_var].ub(); // límite superior de la función objetivo evaluada en la caja actual
		double bigger_diam = thebuffer->top()->box.max_diam(); // dentro de todas las variables (x_i) guarda el diámetro más grande
		double gap_rel_loup = calculate_improved_gap_rel_loup(current_loup, lb_f_obj, ub_f_obj); // gap (resta)
		//int depth = current_cell->depth;

		// Preprocesamos las características
		lb_f_obj = preprocess_feature(lb_f_obj, "lb_f_obj");
		ub_f_obj = preprocess_feature(ub_f_obj, "ub_f_obj");
		bigger_diam = preprocess_feature(bigger_diam, "bigger_diam");
		
		/**********************************************************
		 * FIN Aqui se vienen cambios para comunicar los archivos *
		 **********************************************************/
		
		/**********************************************
		 * AQUI CREACIÓN DE LOS DATOS PARA EL ARCHIVO *
		 **********************************************/
		
		// construimos el JSON
		std::string input_json = "{\"id\": " + std::to_string(call_id++) +
		", \"features\": [";
		input_json += std::to_string(variables) + ", ";
		//input_json += std::to_string(restricciones) + ", ";
		input_json += std::to_string(lb_f_obj) + ", ";
		input_json += std::to_string(ub_f_obj) + ", ";
		input_json += std::to_string(bigger_diam) + ", ";
		input_json += std::to_string(gap_rel_loup);
		//input_json += std::to_string(depth);
		input_json += "]}";
		/**************************************************
		 * FIN AQUI CREACIÓN DE LOS DATOS PARA EL ARCHIVO *
		 **************************************************/

		/*************************************************************
		 * AQUI ENVIAMOS INFORMACIÓN AL MODELO Y RECIBIMOS RESPUESTA *
		 *************************************************************/
			
		python_timer.start();
		std::string response_json = call_python_model(input_json);
		python_timer.stop();
		python_time += python_timer.get_time(); // Simulamos el tiempo de respuesta del modelo
		int model_decision = 0; 
		size_t pos = response_json.find("\"decision\":");
		if (pos != std::string::npos) {
			try {
				model_decision = std::stoi(response_json.substr(pos + 11));
			} catch (const std::exception& e) {
				std::cerr << "Error parsing model decision: " << e.what() << std::endl;
				model_decision = 0; // default to 0 if parsing fails
			}
		}
		
		// Definition of the bisector to be used. 
		// By default, lsmear is used.
		// This can be changed based on the model's decision.
		Bsc* chosen_bsc = &bsc;
		
		switch (model_decision) {
			case 0: chosen_bsc = &bsc; break;
			case 1: chosen_bsc = &bisector_olf; break;
			case 2: chosen_bsc = &bisector_rr; break;
			case 3: chosen_bsc = &bisector_sm; break;
			case 4: chosen_bsc = &bisector_ss; break;
			case 5: chosen_bsc = &bisector_ssr; break;
			default: chosen_bsc = &bsc; break;
		}

		/*****************************************************************
		 * FIN AQUI ENVIAMOS INFORMACIÓN AL MODELO Y RECIBIMOS RESPUESTA *
		 *****************************************************************/
		
		std::ofstream info_file("/home/felipe/Desktop/reuniones_tesis/test/bisectors_probe_nuevas_features.txt", std::ios::app);	
		if (info_file.is_open())
			info_file << model_decision << " bisector inicial " << pid1 << std::endl;

		// cerramos la conexión con el modelo, ya que consultamos en el primer nodo
		//close_model_connection();

		//info_file.close();

		/******************************************
		 * COMIENZA PROCEDIMIENTO DE OPTIMIZACIÓN *
		 ******************************************/
		while (!thebuffer->empty()) {

			loup_changed=false;
			// for double heap , choose randomly the buffer : top  has to be called before pop
			Cell *c = thebuffer->top();
			if (trace >= 2) std::cout << " current box " << c->box << endl;

			try {
				
				// Comprueba si han pasado más de 1800 segundos (30 minutos)
				auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_wall_time).count();
				if (elapsed >= 1800) {
					//std::cout << "Tiempo de ejecución excedido (30 minutos). Terminando optimización." << std::endl;
					status = TIME_OUT;
					// Guardamos el tiempo transcurrido como wall_time
					wall_time = elapsed;
					break;
				}
				
				/**
				 * FUTUREBUFFER ES 0? SI LO ES, LLAMAR AL MODELO, SINO, SIGUE EL MISMO
				 */
				
				if(thebuffer->futurebuffer.size() == 0){
				  	cont++;
				
				// 	/*************************************************************
				// 	 * AQUI ENVIAMOS INFORMACIÓN AL MODELO Y RECIBIMOS RESPUESTA *
				// 	 *************************************************************/
					if(cont % 100 == 0) { // cada 100 iteraciones llamamos al modelo
						variables = n;
						//restricciones = lfd->finder_x_taylor.sys.nb_ctr;
						lb_f_obj = (c->box)[goal_var].lb();
						ub_f_obj = (c->box)[goal_var].ub();
						bigger_diam = c->box.max_diam();
						current_loup = loup;
						gap_rel_loup = calculate_improved_gap_rel_loup(current_loup, lb_f_obj, ub_f_obj);
						// current_cell = c;
						// depth = c->depth;

						std::string input_data_without_preprocessing = "{\"id\": " + std::to_string(call_id) +
							", \"features\": [";
							input_data_without_preprocessing += std::to_string(variables) + ", ";
							//input_data_without_preprocessing += std::to_string(restricciones) + ", ";
							input_data_without_preprocessing += std::to_string(lb_f_obj) + ", ";
							input_data_without_preprocessing += std::to_string(ub_f_obj) + ", ";
							input_data_without_preprocessing += std::to_string(bigger_diam) + ", ";
							input_data_without_preprocessing += std::to_string(gap_rel_loup);
							//input_data_without_preprocessing += std::to_string(depth);
							input_data_without_preprocessing += "]}";
					// 	/**************************************************
					// 	 * FIN AQUI CREACIÓN DE LOS DATOS PARA EL ARCHIVO *
					// 	 **************************************************/

					// 	std::cout << "Input JSON: " << input_data_without_preprocessing << std::endl;
		                
					// 	// cout << "ub - lb: " << ub_f_obj - lb_f_obj << endl;
					// 	// cout << "depth " << depth << endl;

						lb_f_obj = preprocess_feature(lb_f_obj, "lb_f_obj");
						ub_f_obj = preprocess_feature(ub_f_obj, "ub_f_obj");
						bigger_diam = preprocess_feature(bigger_diam, "bigger_diam");

				 	// 	/***************************************
				 	// 	 * FIN CAMBIOS PARA COMUNICAR ARCHIVOS *
				 	// 	 ***************************************/
		
				 	// 	/**********************************************
				 	// 	 * AQUI CREACIÓN DE LOS DATOS PARA EL ARCHIVO *
				 	// 	 **********************************************/
				 	// 	//construimos el JSON
						std::string input_json = "{\"id\": " + std::to_string(call_id++) +
							", \"features\": [";
							input_json += std::to_string(variables) + ", ";
							//input_json += std::to_string(restricciones) + ", ";
							input_json += std::to_string(lb_f_obj) + ", ";
							input_json += std::to_string(ub_f_obj) + ", ";
							input_json += std::to_string(bigger_diam) + ", ";
							input_json += std::to_string(gap_rel_loup);
							//input_json += std::to_string(depth);
							input_json += "]}";

						std::cout << "json: " << input_json << std::endl;
					// 	/**************************************************
					// 	 * FIN AQUI CREACIÓN DE LOS DATOS PARA EL ARCHIVO *
					// 	 **************************************************/
					
						python_timer.start();
						std::string response_json = call_python_model(input_json);
						python_timer.stop();
						python_time += python_timer.get_time(); // Simulamos el tiempo de respuesta del modelo
						int model_decision = 0; 
						size_t pos = response_json.find("\"decision\":");
						if (pos != std::string::npos) {
							try {
								model_decision = std::stoi(response_json.substr(pos + 11));
							} catch (const std::exception& e) {
								std::cerr << "Error parsing model decision: " << e.what() << std::endl;
								model_decision = 0; // default to 0 if parsing fails
							}
						}

					//  	//0: lsmear, 1: lf, 2: rr, 3: sm, 4: ss, 5: ssr
						switch (model_decision) {
							case 0: chosen_bsc = &bsc; break;
							case 1: chosen_bsc = &bisector_olf; break;
							case 2: chosen_bsc = &bisector_rr; break;
							case 3: chosen_bsc = &bisector_sm; break;
							case 4: chosen_bsc = &bisector_ss; break;
							case 5: chosen_bsc = &bisector_ssr; break;
							default: chosen_bsc = &bsc; break;
						}

					// 	if (info_file.is_open())
					// 		info_file << model_decision << " " << pid1 << std::endl;					
					
				 	} //este pertenece al if del contador de feasible divings
				} //este pertenece al if del futurebuffer.size() == 0
				/*****************************************************************
				 * FIN AQUI ENVIAMOS INFORMACIÓN AL MODELO Y RECIBIMOS RESPUESTA *
				 *****************************************************************/

				/******************
				 * METODO ONLINE. *
				 ******************/
				
				pair<Cell*, Cell*> new_cells = chosen_bsc->bisect(*c);

				/**********************
				 * FIN METODO ONLINE. *
				 **********************/
 
				thebuffer->pop();
				delete c; // deletes the cell.

				nb_cells+=2;  // counting the cells handled ( in previous versions nb_cells was the number of cells put into the buffer after being handled)

				handle_cell(*new_cells.first);
				handle_cell(*new_cells.second);

				if (uplo_of_epsboxes == NEG_INFINITY) {
					break;
				}

				if (loup_changed) {
					// In case of a new upper bound (loup_changed == true), all the boxes
					// with a lower bound greater than (loup - goal_prec) are removed and deleted.
					// Note: if contraction was before bisection, we could have the problem
					// that the current cell is removed by contractHeap. See comments in
					// older version of the code (before revision 284).

					double ymax=compute_ymax();

					thebuffer->contract(ymax);

					//cout << " now buffer is contracted and min=" << buffer.minimum() << endl;

					// TODO: check if happens. What is the return code in this case?
					if (ymax <= NEG_INFINITY) {
						if (trace) std::cout << " infinite value for the minimum " << endl;
						break;
					}
				}
				update_uplo();

				if (!anticipated_upper_bounding) // useless to check precision on objective if 'true'
					if (get_obj_rel_prec()<rel_eps_f || get_obj_abs_prec()<abs_eps_f)
						break;

				if (timeout>0) timer.check(timeout); // TODO: not reentrant, JN: done
				time = timer.get_time();

			} catch (NoBisectableVariableException& ) {
				update_uplo_of_epsboxes((c->box)[goal_var].lb());
				thebuffer->pop();
				delete c; // deletes the cell.
				update_uplo(); // the heap has changed -> recalculate the uplo (eg: if not in best-first search)

			}
		}
		/*******************************************
		 *** ACABA PROCEDIMIENTO DE OPTIMIZACIÓN ***
		 *******************************************/
		
		timer.stop();
		time = timer.get_time();
		
		//calculamos el tiempo de ejecución real/total
		auto end_wall_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration<double>(end_wall_time - start_wall_time);
		wall_time = duration.count();
		
		//info_file.close();
		//std::cout << "# fd: " << cont << endl;

		// Cerramos la conexión con el modelo
		close_model_connection();

		// No solution found and optimization stopped with empty buffer
		// before the required precision is reached => means infeasible problem
	 	if (uplo_of_epsboxes == NEG_INFINITY)
	 		status = UNBOUNDED_OBJ;
	 	else if (uplo_of_epsboxes == POS_INFINITY && (loup==POS_INFINITY || (loup==initial_loup && abs_eps_f==0 && rel_eps_f==0)))
	 		status = INFEASIBLE;
	 	else if (loup==initial_loup)
	 		status = NO_FEASIBLE_FOUND;
	 	else if (get_obj_rel_prec()>rel_eps_f && get_obj_abs_prec()>abs_eps_f)
	 		status = UNREACHED_PREC;
	 	else
	 		status = SUCCESS;
	} catch (TimeOutException& ) {
		status = TIME_OUT;
	}

	/* TODO: cannot retrieve variable names here. */
	for (int i=0; i<(extended_COV ? n+1 : n); i++)
		cov->data->_optim_var_names.push_back(string(""));

	cov->data->_optim_optimizer_status = (unsigned int) status;
	cov->data->_optim_uplo = uplo;
	cov->data->_optim_uplo_of_epsboxes = uplo_of_epsboxes;
	cov->data->_optim_loup = loup;

	cov->data->_optim_time += time;
	cov->data->_optim_nb_cells += nb_cells;
	cov->data->_optim_loup_point = loup_point;

	// for conversion between original/extended boxes
	IntervalVector tmp(extended_COV ? n+1 : n);

	// by convention, the first box has to be the loup-point.
	if (extended_COV) {
		write_ext_box(loup_point, tmp);
		tmp[goal_var] = Interval(uplo,loup);
		cov->add(tmp);
	}
	else {
		cov->add(loup_point);
	}

	while (!buffer.empty()) {
		Cell* cell=buffer.top();
		if (extended_COV)
			cov->add(cell->box);
		else {
			read_ext_box(cell->box,tmp);
			cov->add(tmp);
		}
		delete buffer.pop();
	}

	return status;
}

namespace {
const char* green() {
#ifndef _WIN32
	return "\033[32m";
#else
	return "";
#endif
}

const char* red(){
#ifndef _WIN32
	return "\033[31m";
#else
	return "";
#endif
}

const char* white() {
#ifndef _WIN32
	return "\033[0m";
#else
	return "";
#endif
}

}

void Optimizer::report() {

	double ibex_time = time - python_time;
	
	/**
	 * obtenemos el pid
	 */

	pid_t pid = getpid();

	// if (!cov || !buffer.empty()) { // not started
	// 	cout << " not started." << endl;
	// 	return;
	// }

	// switch(status) {
	// case SUCCESS:
	// 	cout << green() << " optimization successful!" << endl;
	// 	break;
	// case INFEASIBLE:
	// 	cout << red() << " infeasible problem" << endl;
	// 	break;
	// case NO_FEASIBLE_FOUND:
	// 	cout << red() << " no feasible point found (the problem may be infeasible)" << endl;
	// 	break;
	// case UNBOUNDED_OBJ:
	// 	cout << red() << " possibly unbounded objective (f*=-oo)" << endl;
	// 	break;
	// case TIME_OUT:
	// 	cout << red() << " time limit " << timeout << "s. reached " << endl;
	// 	break;
	// case UNREACHED_PREC:
	// 	cout << red() << " unreached precision" << endl;
	// 	break;
	// }
	// cout << white() <<  endl;

	// // No solution found and optimization stopped with empty buffer
	// // before the required precision is reached => means infeasible problem
	// if (status==INFEASIBLE) {
	// 	cout << " infeasible problem " << endl;
	// } else {
	// 	cout << " f* in\t[" << uplo << "," << loup << "]" << endl;
	// 	cout << "\t(best bound)" << endl << endl;

	// 	if (loup==initial_loup)
	// 		cout << " x* =\t--\n\t(no feasible point found)" << endl;
	// 	else {
	// 		if (loup_finder.rigorous())
	// 			cout << " x* in\t" << loup_point << endl;
	// 		else
	// 			cout << " x* =\t" << loup_point.lb() << endl;
	// 		cout << "\t(best feasible point)" << endl;
	// 	}
	// 	cout << endl;
	// 	double rel_prec=get_obj_rel_prec();
	// 	double abs_prec=get_obj_abs_prec();

	// 	cout << " relative precision on f*:\t" << rel_prec;
	// 	if (rel_prec <= rel_eps_f)
	// 		cout << green() << " [passed] " << white();
	// 	cout << endl;

	// 	cout << " absolute precision on f*:\t" << abs_prec;
	// 	if (abs_prec <= abs_eps_f)
	// 		cout << green() << " [passed] " << white();
	// 	cout << endl;
	// }

	// cout << " cpu time used:\t\t\t" << time << "s";
	// if (cov->time()!=time)
	// cout << " [total=" << cov->time() << "]";
	// cout << endl;

	// // ==========================================================
    // // INICIO DEL NUEVO REPORTE DE TIEMPO
    // // ==========================================================
    // cout << " tiempo en llamadas a Python:\t" << python_time << "s";
    // if (time > 0) {
    //     cout << " (" << fixed << setprecision(2) << (python_time / time) * 100.0 << "%)";
    // }
    // cout << endl;
    // ==========================================================
    // FIN DEL NUEVO REPORTE DE TIEMPO
    // ==========================================================
	
	// cout << " number of cells:\t\t" << nb_cells;
	// if (cov->nb_cells()!=nb_cells)
	// 	cout << " [total=" << cov->nb_cells() << "]";
	// cout << endl << endl;
	
	// if (statistics) 
	// 	cout << "  ===== Statistics ====" << endl << endl << *statistics << endl;

	// cout << nb_cells << " " << time << endl;
	
	//cout << nb_cells << " " << time << " " << python_time << " " << ibex_time << " " << wall_time << " " << endl;
	 
	std::cout << nb_cells << " " << wall_time << endl;

	// registramos el tiempo real en un archivo
	// std::ofstream time_file("/home/felipe/Desktop/reuniones_tesis/test/wall_time_probe.txt", std::ios::app);
	// if (time_file.is_open()) {
	// 	time_file << wall_time << " " << pid << std::endl;
	// 	time_file.close();
	// }
	
}



} // end namespace ibex
