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
#include <iomanip>
#include <queue>
#include <string>
#include <iostream>
#include <fstream>

#include <limits> // Necesario para std::numeric_limits
#include <cmath>  // Necesario para std::abs


// ============================================================================
//         COMUNICACIÓN CON MODELO MEDIANTE PIPES
// ============================================================================

// declaramos aqui las rutas de las tuberias de comunicacion
const std::string REQUEST_PIPE_PATH = "/tmp/ibex_request_pipe";
const std::string RESPONSE_PIPE_PATH = "/tmp/ibex_response_pipe";

// realizamos la nueva funcion de comunicacion entre c++ y python
std::string call_python_model(const std::string& input_json) {
    // abrimos la tuberia para escribir (solicitud)
    // el constructor de ofstream va a abrirla
    std::ofstream req_pipe(REQUEST_PIPE_PATH);
    if (!req_pipe.is_open()) {
        std::cerr << "Error: Could not open request pipe." << std::endl;
        // return "{\"error\": \"Could not open request pipe\"}";
    }
    
    // escribimos el JSON seguido de un salto de linea y cerramos
    // el flush es importante para asegurar que se envie de manera inmediata
    // cerrar la tuberia es la señal para que lector de python sepa que terminanos
    // de escribir
    req_pipe << input_json << std::endl;
    req_pipe.close();

    // abrimos la tuberia para leer (respuesta)
    // esta llamada se bloquea hasta que el servidor escriba la respuesta
    std::ifstream res_pipe(RESPONSE_PIPE_PATH);
    if (!res_pipe.is_open()) {
        std::cerr << "Error: Could not open response pipe." << std::endl;
       // return "{\"error\": \"Could not open response pipe\"}";
    }
    
    // leemos la respuesta (una sola linea)
    std::string response_str;
    if (!std::getline(res_pipe, response_str)) {
        std::cerr << "Error: Could not read from response pipe." << std::endl;
        // return "{\"error\": \"Failed to read from response pipe\"}";
    }
    
    res_pipe.close();
    
    return response_str;
}

// ============================================================================
//                      FIN DE LA COMUNICACIÓN CON EL MODELO
// ============================================================================

// ============================================================================
//  AQUÍ SE REALIZA EL PREPROCESAMIENTO DE LOS DATOS ANTES DE LLAMAR AL MODELO
// ============================================================================

double const UMBRAL_CERO = 1e-12; 
double const VALOR_MAXIMO_JSON = 1e30; // Valor seguro para JSON

// Función ÚNICA y GENÉRICA de preprocesamiento
// Mantenemos el argumento 'feature_name' para que tu llamada en build_feature_json no falle,
// aunque internamente ya no necesitemos usarlo para switch/case.
double preprocess_feature(double val, const std::string& feature_name) {
    
    // 1. Manejar NaN: Reemplazar por 0 es seguro
    if (std::isnan(val)) {
        return 0.0; 
    }
    
    // 2. Manejar Infinitos:
    // JSON no soporta "Infinity". Lo reemplazamos por el máximo double seguro.
    if (std::isinf(val)) {
        return (val > 0) ? VALOR_MAXIMO_JSON : -VALOR_MAXIMO_JSON;
    }

    // 3. Limpieza de ruido numérico
    if (std::abs(val) < UMBRAL_CERO) {
        val = 0.0;
    }

    // --- SIN HARD CLIPPING ---
    // Dejamos pasar el valor real (ej: 1e15) para que Python aplique SymLog.
    
    return val;
}

// ================================================================================
//  FIN AQUÍ SE REALIZA EL PREPROCESAMIENTO DE LOS DATOS ANTES DE LLAMAR AL MODELO
// ================================================================================

using namespace std;

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

	if (cov) delete cov;
	cov = new CovOptimData(extended_COV? n+1 : n, extended_COV);
	cov->data->_optim_time = data.time();
	cov->data->_optim_nb_cells = data.nb_cells();
}

std::string build_feature_json(long call_id, const IntervalVector& box, const System& sys, int goal_var) {
    int variables = sys.nb_var;
    int restricciones = sys.nb_ctr;
    
    // Obtener valores crudos
    double lb_f_obj = box[goal_var].lb();
    double ub_f_obj = box[goal_var].ub();
    double search_space = 0.0;
    for (int i = 0; i < sys.nb_var; i++) {
        search_space += box[i].diam();
    }
    double bigger_diam = box.max_diam();
    double lower_diam = box.min_diam();

    // Sanitizar (Solo quitar NaNs/Infs, sin recortar magnitud)
    lb_f_obj = preprocess_feature(lb_f_obj, "lb_f_obj");
    ub_f_obj = preprocess_feature(ub_f_obj, "ub_f_obj");
    search_space = preprocess_feature(search_space, "search_space"); // Asegúrate de procesar este también
    bigger_diam = preprocess_feature(bigger_diam, "bigger_diam");
    lower_diam = preprocess_feature(lower_diam, "lower_diam");

    // Construir JSON
    // NOTA: Asegúrate de que el orden coincida con los índices de Python
    // 0:variables, 1:restricciones, 2:lb, 3:ub, 4:search_space, 5:bigger, 6:lower
    std::string json = "{\"id\": " + std::to_string(call_id) +
                       ", \"features\": [";
    json += std::to_string(variables) + ", ";
    json += std::to_string(restricciones) + ", ";
    json += std::to_string(lb_f_obj) + ", ";
    json += std::to_string(ub_f_obj) + ", ";
    json += std::to_string(search_space) + ", ";
    json += std::to_string(bigger_diam) + ", ";
    json += std::to_string(lower_diam) + "]}";

    return json;
}

Optimizer::Status Optimizer::optimize() {
	Timer timer;
	timer.start();
	
	update_uplo();
    pid_t pid_t = getpid();

    std::ofstream info_file("/home/felipe/Desktop/bisectores_modelo_feasible_diving_pipes/validacion_200sim_online_symlog.txt", std::ios::app);


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
        long call_id = 0;
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

        Bsc* bisectors[] = {&bsc, &bisector_olf, &bisector_rr, &bisector_sm, &bisector_ss, &bisector_ssr};
		Bsc *chosen_bsc = bisectors[0];
        int model_decision = -1;
        bool root_node = true;

        while (!thebuffer->empty()) {

			loup_changed=false;
			// for double heap , choose randomly the buffer : top  has to be called before pop
			Cell *c = thebuffer->top();
			if (trace >= 2) std::cout << " current box " << c->box << endl;

            IntervalVector& box = c->box;

			try {
				
				if(thebuffer->futurebuffer.size() == 0){
				  	cont++;
                }
				
                // cada 500 feasible divings llamamos al modelo (o en el nodo raíz)
                if (cont % 200 == 0 || root_node) { //
                    root_node = false; 
                    std::cout << "id: " << call_id << std::endl;
                    std::string input_json = build_feature_json(call_id++, box, system, goal_var);
                    std::string response_json = call_python_model(input_json);
                    //std::cout << "data: " << input_json << std::endl;
                    //model_decision = 0; 
                    size_t pos = response_json.find("\"decision\":");
                    if (pos != std::string::npos) {
                        try {
                            model_decision = std::stoi(response_json.substr(pos + 11));
                            // std::cout << "response: " << response_json << std::endl;
                        } catch (const std::exception& e) {
                            std::cerr << "Error parsing model decision: " << e.what() << std::endl;
                            model_decision = 0; // default to 0 if parsing fails
                        }
                    }

                    chosen_bsc = bisectors[model_decision];

                    if (info_file.is_open()) {
                        info_file << model_decision << " " << pid_t << std::endl;
                    }
                }
				
				pair<Cell*, Cell*> new_cells = chosen_bsc->bisect(*c);
 
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

    std::ofstream info_file("/home/felipe/Desktop/bisectores_modelo_feasible_diving_pipes/validacion_200sim_online_symlog.txt", std::ios::app);
	pid_t pid_t = getpid();

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
	 
	std::cout << nb_cells << " " << time << endl;
    info_file << nb_cells << " " << time << " " << pid_t << std::endl;
}



} // end namespace ibex
