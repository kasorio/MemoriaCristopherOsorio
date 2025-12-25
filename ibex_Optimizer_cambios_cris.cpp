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
#include <stdlib.h>
#include "ibex_CellBeamSearch.h"
#include "ibex_LoupFinderDefault.h"

#include "ibex_SmearFunction.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_OptimLargestFirst.h"
#include "ibex_System.h"

#include <float.h>
#include <vector>
#include <queue>
#include <stdlib.h>
#include <iomanip>

#include <filesystem> 

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
										time(0), nb_cells(0), cov(NULL) {

	if (trace) cout.precision(12);
	
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
			cout << "                    ";
			cout << "\033[32m loup= " << loup << "\033[0m" << endl;
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
			cout << " loup = " << loup << " new_uplo=" << new_uplo <<  " uplo_of_epsboxes=" << uplo_of_epsboxes << endl;
			ibex_error("optimizer: new_uplo>loup (please report bug)");
		}
		if (new_uplo < uplo) {
			cout << "uplo= " << uplo << " new_uplo=" << new_uplo << endl;
			ibex_error("optimizer: new_uplo<uplo (please report bug)");
		}

		// uplo <- max(uplo, min(new_uplo, uplo_of_epsboxes))
		if (new_uplo < uplo_of_epsboxes) {
			if (new_uplo > uplo) {
				uplo = new_uplo;

				if (trace)
					cout << "\033[33m uplo= " << uplo << "\033[0m" << endl;
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
			cout << " unprocessable tiny box: now uplo<=" << setprecision(12) <<  uplo_of_epsboxes << " uplo=" << uplo << endl;
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

	loup=obj_init_bound; // loup => lower upperbound, uplo => upper lowerbound


	/***************************
	 ** INICIO MODIFICACIONES **
	 ***************************/

	// double search_space = 1;
	// double bigger_diam = NEG_INFINITY;
	// double lower_diam = POS_INFINITY;
	
	IntervalVector aux(init_box.size());  //crea una variable auxiliar del tamaño de la cantidad de variables

	for (int i = 0; i < init_box.size(); i++) {
		aux[i] = init_box[i]; //asigna el valor de init_box a la variable auxiliar

		//+10% 
		//aux[i] = Interval((aux[i].lb() - aux[i].diam()/20), (aux[i].ub() + aux[i].diam()/20)); // esta aumenta un 10%

		//+25% 
		//aux[i] = Interval((aux[i].lb() - aux[i].diam()/8), (aux[i].ub() + aux[i].diam()/8)); // esta aumenta un 25%

		//+50% 
		//aux[i] = Interval((aux[i].lb() - aux[i].diam()/4), (aux[i].ub() + aux[i].diam()/4)); // esta aumenta un 50%
	}
	
	// Just to initialize the "loup" for the buffer
	// TODO: replace with a set_loup function
	
	buffer.contract(loup);

	uplo=NEG_INFINITY;
	uplo_of_epsboxes=POS_INFINITY;

	nb_cells=0;

	buffer.flush();

	Cell* root=new Cell(IntervalVector(n+1));

	write_ext_box(aux, root->box);

	// add data required by the bisector
	bsc.add_property(aux, root->prop);

	// add data required by the contractor
	ctc.add_property(aux, root->prop);

	// add data required by the buffer
	buffer.add_property(aux, root->prop);

	// add data required by the loup finder
	loup_finder.add_property(aux, root->prop);

	//cout << "**** Properties ****\n" << root->prop << endl;

	loup_changed=false;
	initial_loup=obj_init_bound;

	loup_point = aux; //.set_empty();
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

/**
 * CALCULO DE LA FEATURE GAP REL LOUP
 * EL FEATURE ES MÁS ROBUSTO QUE EL ANTERIOR, POR LO QUE DEBIESE FUNCIONAR MEJOR
 * \brief Calculate a more robust relative gap feature based on the loup.
 */

// double calculate_improved_gap_rel_loup(double current_loup, double lb_f_obj, double ub_f_obj) {
//     double epsilon = 1e-12; // Valor más pequeño para mayor precisión
//     double gap_rel_loup;
    
//     if (!std::isfinite(current_loup) || current_loup == POS_INFINITY) {
//         // Caso 1: No hay loup válido
//         // El gap representa la incertidumbre relativa de la caja actual
//         double numerator = ub_f_obj - lb_f_obj;
//         double denominator = std::max(fabs(ub_f_obj), fabs(lb_f_obj));
        
//         if (denominator < epsilon) {
//             gap_rel_loup = 1.0; // Máxima incertidumbre cuando ambos límites son ~0
//         } else {
//             gap_rel_loup = numerator / (denominator + epsilon);
//         }
        
//         // Normalizar a [0,1] - mayor valor = menos prometedora
//         gap_rel_loup = std::min(1.0, std::max(0.0, gap_rel_loup));
        
//     } else {
//         // Caso 2: Hay loup válido
//         // Calculamos qué tan cerca está el lower bound de la mejor solución conocida
//         double numerator = current_loup - lb_f_obj;
        
//         // Si lb_f_obj > current_loup, la caja es muy prometedora (gap negativo)
//         if (numerator < 0) {
//             gap_rel_loup = 0.0; // Muy prometedora
//         } else {
//             // Normalizamos por la magnitud del loup
//             double denominator = fabs(current_loup) + epsilon;
//             gap_rel_loup = numerator / denominator;
            
//             // Limitamos a [0,1] para evitar valores extremos
//             gap_rel_loup = std::min(1.0, gap_rel_loup);
//         }
//     }
    
//     return gap_rel_loup;
// }

/**
 * Feature adicional: qué tan prometedora es la caja comparada con otras en el buffer
 * Esta feature es de prueba por ahora, solo para ver como funciona
 */
// double calculate_relative_promise(double lb_f_obj, double buffer_minimum) {
//     double epsilon = 1e-12;
    
//     if (!std::isfinite(buffer_minimum)) {
//         return 0.5; // Valor neutro si no hay referencia
//     }
    
//     double diff = lb_f_obj - buffer_minimum;
    
//     if (fabs(diff) < epsilon) {
//         return 1.0; // Es la más prometedora del buffer
//     }
    
//     // Normalizamos: valores más cercanos a buffer_minimum son más prometedores
//     double scale = std::max(fabs(buffer_minimum), 1.0);
//     return std::exp(-fabs(diff) / scale); // Decay exponencial
// }

double calculate_bounds_ratios(double lb, double ub, double epsilon, int id) {
    static double initial_lb, initial_ub, initial_obj_range;
    static bool initial_bounds_set = false;
    
    // Si es la primera llamada (id == 1), guardar valores iniciales
    if (id == 1) {
        initial_lb = lb;
        initial_ub = ub;
        initial_obj_range = initial_ub - initial_lb;
        initial_bounds_set = true;
    }
    
    // Verificar que se hayan establecido los valores iniciales
    if (!initial_bounds_set) {
        return -5.0; // Error: no se han establecido valores iniciales
    }
    
    // Calcular rango actual
    double current_obj_range = ub - lb;
    
    // Edge cases para bounds iniciales
    bool initial_bounds_infinite = (!std::isfinite(initial_lb) || !std::isfinite(initial_ub));
    bool initial_bounds_invalid = (!std::isfinite(initial_obj_range) || std::abs(initial_obj_range) < epsilon);
    
    // Edge cases para bounds actuales
    bool current_bounds_invalid = (!std::isfinite(current_obj_range) || std::abs(current_obj_range) < epsilon);
    
    // Caso 1: Current bounds inválidos
    if (current_bounds_invalid) {
        return -5.0; // Valor sentinel
    }
    
    // Caso 2: Bounds iniciales problemáticos - usar solo rango actual
    if (initial_bounds_infinite || initial_bounds_invalid) {
        double log_current_range = std::log10(current_obj_range);
        return std::max(-5.0, std::min(5.0, log_current_range));
    }
    
    // Caso 3: Normal - calcular progreso relativo
    double progress_ratio = current_obj_range / initial_obj_range;
    return std::max(0.0, std::min(1.0, progress_ratio));
}

double calculate_diameter_shape(double bigger_diam, double lower_diam, double epsilon) {
    // Edge cases
    if (lower_diam < epsilon || !std::isfinite(lower_diam)) return 0.0;
    if (bigger_diam < epsilon || !std::isfinite(bigger_diam)) return 0.0;
    
    double diam_ratio = bigger_diam / lower_diam;
    if (diam_ratio <= 1.0 + epsilon) return 0.0; // Box cuadrado
    
    double log_diam_ratio = std::log10(diam_ratio);
    return std::max(-3.0, std::min(3.0, log_diam_ratio));
}

double calculate_diameter_progress(double current_diam, double epsilon, int id, bool is_bigger_diam) {
    static double initial_bigger_diam, initial_lower_diam;
    static bool initial_set = false;
    
    if (id == 1) {
        if (is_bigger_diam) {
            initial_bigger_diam = current_diam;
        } else {
            initial_lower_diam = current_diam;
        }
        initial_set = true;
        return 1.0; // Valor inicial
    }
   
    if (!initial_set) return 0.5;
    
    double initial_value = is_bigger_diam ? initial_bigger_diam : initial_lower_diam;
    if (initial_value < epsilon) return 0.5;
    
    double progress = current_diam / initial_value;
    return std::max(0.0, std::min(1.0, progress));
} 

Optimizer::Status Optimizer::optimize() {
	Timer timer;
	timer.start();

	update_uplo();

	try {

		/*****************************************
		 ** DECLARACIÓN DE BUFFER Y LOUP FINDER **
		 *****************************************/

		/**
		 * \note dynamic_cast convierte de manera "forzosa" el puntero de la clase base a un puntero de la clase derivada.
		 */ 

		CellBeamSearch * thebuffer = dynamic_cast<CellBeamSearch*>(&buffer);
		LoupFinderDefault * lfd = dynamic_cast<LoupFinderDefault*>(&loup_finder);

		//vector of cells
		queue<Cell*> aux; 

		//variable auxiliar para guardar el valor de la caja "inicial"
		IntervalVector aux_box = thebuffer->top()->box;
		double prec = 1e-7;

		/***********************************
		 ** DECLARACIÓN DE LOS BISECTORES **
		 ***********************************/
		OptimLargestFirst bisector_olf(goal_var, true, prec, 0.5);
		RoundRobin bisector_rr(prec, 0.5);

		// obtención del sistema para bisectores smear
		System system = lfd->finder_x_taylor.sys;
		SmearMax bisector_sm(system,prec);
		SmearSum bisector_ss(system,prec);
		SmearSumRelative bisector_ssr(system,prec);
		
		// upper lowerbound
		double aux_uplo = uplo;
		// lower upperbound
		double aux_loup = loup;
		// región interior donde se saca el uplo y loup
		IntervalVector inner=aux_box;
		
		//		while(!buffer.empty())
		// simulaciones son la cantidad de datos que se generarán
		int num_sim = 1000;
		double epsilon = 1e-9;

		// se ingresa al vector la caja inicial del buffer (primer nodo a tratar)
		aux.push(thebuffer->top());
		
		// se obtiene la celda actual para el análisis
		// aqui calculamos los bounds iniciales para posteriormente calcular una proporción
		double lb_f_obj_inicial = lfd->finder_x_taylor.sys.goal->eval(aux.front()->box).lb();
		double ub_f_obj_inicial = lfd->finder_x_taylor.sys.goal->eval(aux.front()->box).ub();
		double bigger_diam_inicial = aux.front()->box.max_diam();
		double lower_diam_inicial = aux.front()->box.min_diam();		

		std::cout << "lb_f_obj_inicial: " << lb_f_obj_inicial << std::endl;
		std::cout << "ub_f_obj_inicial: " << ub_f_obj_inicial << std::endl;
		std::cout << "bigger_diam_inicial: " << bigger_diam_inicial << std::endl;
		std::cout << "lower_diam_inicial: " << lower_diam_inicial << std::endl;
		// estas variables son las que utilizaremos para calcular los límites dentro de las simulaciones.
		double lb_f_obj;
		double ub_f_obj;
		double bigger_diam;
		double lower_diam;
		
		// Variables finales que se utilizarán como input en la red neuronal
		int variables = n;
		BitSet active;
		double ratio_bounds;
		double ratio_bigger_diam;
		double ratio_lower_diam;
		double box_shape;
		
		// Obtener información del benchmark desde variable de entorno
		char* benchmark_env = getenv("CURRENT_BENCHMARK");
		string benchmark_name = "unknown";
		if (benchmark_env) {
			benchmark_name = string(benchmark_env);
		}

		// Crear nombres de archivo basados en el benchmark
		string input_filename = "/home/cristopher/Escritorio/dataset/input/optimizer_input_" + benchmark_name + ".txt";
		string output_filename = "/home/cristopher/Escritorio/dataset/output/optimizer_output_" + benchmark_name + ".txt";

		// Crear directorios si no existen
		std::system("mkdir -p /home/cristopher/Escritorio/dataset/input");
		std::system("mkdir -p /home/cristopher/Escritorio/dataset/output");

		std::ofstream InputFile(input_filename, std::ios::app);
		std::ofstream OutputFile(output_filename, std::ios::app);
		if (!InputFile.is_open()) {
			cerr << "No se pudo abrir el archivo. Comprueba la ruta y permisos." << endl;
			exit(1);
		}

		if (!OutputFile.is_open()) {
			cerr << "No se pudo abrir el archivo. Comprueba la ruta y permisos." << endl;
			exit(1);
		}

		// se ingresa la cantidad de simulaciones que se realizarán
		for (int k = 0 ; k < num_sim ; k++){

            // si el tamaño del vector es 0, se sale del ciclo
            if(aux.size() == 0) break;

            // --- Obtenemos la celda actual (nodo) ---
            Cell* current_cell = aux.front();

            // Feature 1: variables
            int feat_variables = n;

            // Feature 2: restricciones
            BitSet active = lfd->finder_x_taylor.sys.active_ctrs(current_cell->box);
            int feat_restricciones = active.size();

            // Feature 3: profundidad
            // 'current_cell' es un 'CellBeamSearch' que SÍ tiene 'depth'
            int feat_profundidad = current_cell->depth;

            // Feature 4: cota_superior (Loup)
            // Usamos la cota global 'loup' (copiada en 'aux_loup' al inicio)
            double feat_cota_superior = aux_loup; 

            // Feature 5: cota_inferior (Uplo)
            // Usamos la cota global 'uplo' (copiada en 'aux_uplo' al inicio)
            double feat_cota_inferior = aux_uplo;

            // Features 6 y 7: diametro_grande y diametro_pequeno
            // ¡IMPORTANTE! Deben ser los diámetros de las variables 'x',
            // no de la caja extendida 'c->box' (que incluye 'y').
            // Usamos la función 'read_ext_box' para obtener solo las 'x'.
            
            IntervalVector x_box(n); // Un vector para almacenar solo las variables 'x'
            read_ext_box(current_cell->box, x_box); // Llenamos x_box
            
            double feat_diametro_grande = x_box.max_diam();
            double feat_diametro_pequeno = x_box.min_diam();

            // --- ESCRITURA EN ARCHIVO ---
            // Reemplazamos las líneas de escritura anteriores por estas:

            InputFile << "variables: " << feat_variables << endl;
            InputFile << "restricciones: " << feat_restricciones << endl;
            InputFile << "profundidad: " << feat_profundidad << endl;
            InputFile << "cota_superior: " << feat_cota_superior << endl;
            InputFile << "cota_inferior: " << feat_cota_inferior << endl;
            InputFile << "diametro_grande: " << feat_diametro_grande << endl;
            InputFile << "diametro_pequeno: " << feat_diametro_pequeno << endl;
            InputFile << "id: " << k+1 << endl << endl; // Separador para la siguiente entrada


            /***********************************************************************
             * Se escriben en el archivo los datos para el input de la red neuronal *
             ***********************************************************************/
            
            // NOTA: He eliminado las líneas que escribían:
            // ratio_bounds, ratio_bigger_diam, ratio_lower_diam, box_shape
            // Si también quieres conservarlas, solo añádelas aquí.

            // El resto del bucle (la simulación de heurísticas) permanece igual...

			/***********************************************************************
			 * Se escriben en el archivo los datos para el input de la red neuronal *
			 ***********************************************************************/
			
			// para cada técnica (lsmear, rr, olf, smearmax, smearsum, smearsumrel)
			for (int i = 0 ; i < 6 ; i++){

				// se asignan los nuevos valores frontera
				uplo = aux_uplo;
				loup = aux_loup;

				// aun no se han bisectado nodos
				nb_cells = 0;

				// si no hay nodos (es decir, es el primero), se le asigna uno desde el vector auxiliar
				if(thebuffer->size() == 0){
					//cout << aux.size() << endl; exit(1);
					thebuffer->push(aux.front());
				}

				bool first_iteration = true;
				
				// variable para saber si la técnica se registró o no
				bool technique_registered = false;
				
				// mientras queden nodos por revisar
				while (!thebuffer->empty()) {

					// el loup no ha cambiado
					loup_changed=false;

					// for double heap , choose randomly the buffer : top  has to be called before pop
					// celda "padre" (la que se está revisando)
					Cell *c = thebuffer->top();

					if (trace >= 2) std::cout << " current box " << c->box << endl;

					try {
						// se crea un par de celdas que corresponden a los nodos "hijos" con cada técnica
						pair<Cell*,Cell*> new_cells;
						// new_cells = bsc.bisect(*c); 

						// comienza la bisección dependiendo cada técnica
						if (i == 0) //lsmear
							new_cells=bsc.bisect(*c);
						if (i == 1) //lf
							new_cells=bisector_olf.bisect(*c);
						if (i == 2) //rr
							new_cells=bisector_rr.bisect(*c);
						if (i == 3) //smearmax
							new_cells=bisector_sm.bisect(*c);
						if (i == 4) //smearsum
							new_cells=bisector_ss.bisect(*c);
						if (i == 5) //smearsumrel
							new_cells=bisector_ssr.bisect(*c);

						// se elimina el nodo "padre" de la lista de nodos por revisar
						thebuffer->pop();

	// 					// this is part of the modification 
						if (first_iteration){  //to save the reference to the first node
							first_iteration = false;
						}
						
						else {
							delete c; // se elimina el nodo "padre" de la memoria
						} 

						nb_cells+=2;  // counting the cells handled ( in previous versions nb_cells was the number of cells put into the buffer after being handled)
						
						// se manejan las celdas hijas
						handle_cell(*new_cells.first);
						handle_cell(*new_cells.second);

						// se revisa si ya no hay más nodos por revisar
						// thebuffer es el buffer que contiene los nodos próximos a visitar
						// futurebuffer es el buffer que se hace en el feasiblediving
						// si entra aquí, va a buscar la siguiente caja a bisectar

						// si futurebuffer es 0, es porque se acabó la busqueda en profundidad
						// se va a cambiar el nodo, y va a comenzar otra busqueda


						if(thebuffer->get_futurebuffer().size() == 0){ //deadend has arrived
							if (i == 0){
								OutputFile << "LSMEAR ";
								// std::cout << "LSMEAR ";
							}
							else if (i == 1){
								OutputFile <<"LF ";
								// std::cout <<"LF ";
							}
							else if (i == 2){
								OutputFile << "RR ";
								// std::cout << "RR ";
							}
							else if (i == 3){
								OutputFile << "SM ";
								// std::cout << "SM ";
							}
							else if (i == 4){
								OutputFile << "SS ";
								// std::cout << "SS ";
							}
							else if (i == 5){
								OutputFile << "SSR ";
								// std::cout << "SSR ";
							}
							OutputFile << nb_cells << endl;
							technique_registered = true;
							
							// std::cout << nb_cells << endl;

						// 	// REVISAR Y PREGUNTAR AL PROFE QUE HACE ESTO
						// 	// finalmente esto guarda los nodos de la lsmear y elimina los que generan las técnicas (serían nodos thrashing)
							int auxaux=thebuffer->size();
							for (int tt = 0 ; tt < auxaux ; tt++){
								if(i==0){ //se copia el buffer para poder usarlo en las siguientes fases
									Cell *c = thebuffer->top();
									aux.push(c);
									thebuffer->pop();
								} // este if corresponde a los nodos que generan las técnicas (rr, lf, sm, ss, ssr)
								else{
									Cell *c = thebuffer->top();
									thebuffer->pop();
									delete c; // deletes the cell.
								} // esta parte corresponde a los nodos que genera lsmear

						    }
	// //						if(i==0){
	// //
	// //							for (int tt = 0 ; tt < auxaux ; tt++){
	// //								Cell *c = thebuffer->CellHeap::top();
	// //								aux.push(c);
	// //								thebuffer->CellHeap::pop();
	// //							}
	// //						}
	// //						thebuffer->CellHeap::flush();
						} 
							//else {
							// 	cout << " still " << thebuffer->size() << " cells to process" << endl;
							// }
		//					else{
		//						total_nodes = total_nodes+thebuffer->get_futurebuffer().size();
		//					}

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

					}

					catch (NoBisectableVariableException& ) {
						update_uplo_of_epsboxes((c->box)[goal_var].lb());
						thebuffer->pop();
						if(nb_cells!=0)
							delete c; // deletes the cell.
						update_uplo(); // the heap has changed -> recalculate the uplo (eg: if not in best-first search)
					}

				}

				// en caso de no registrarse la técnica en el if del deadend
				// se registra aquí
				// puede ser que estos registros se hagan por término de precisión, factibilidad, etc
				if (!technique_registered){
					if (i == 0){
						OutputFile << "LSMEAR ";
						//std::cout << "LSMEAR no fd ";
						// std::cout << "LSMEAR ";
					}
					else if (i == 1){
						OutputFile <<"LF ";
						//std::cout << "LF no fd ";
						// std::cout <<"LF ";
					}
					else if (i == 2){
						OutputFile << "RR ";
						//std::cout << "RR no fd ";
						// std::cout << "RR ";
					}
					else if (i == 3){
						OutputFile << "SM ";
						//std::cout << "SM no fd ";
					}
					else if (i == 4){
						OutputFile << "SS ";
						//std::cout << "SS no fd ";
					}
					else if (i == 5){
						OutputFile << "SSR ";
						//std::cout << "SSR no fd ";
					}
					OutputFile << nb_cells << endl;
					//std::cout << nb_cells << " " << k+1 << endl;
					// std::cout << nb_cells << endl;
				}
			}
			OutputFile << "id: " << k+1 << endl << endl;
			//std::cout << "id: " << k+1 << endl << endl;
			Cell *c = aux.front();
			aux.pop();
			delete c; // deletes the cell.
		}

		InputFile << "------------------------------------------------"  << endl;
		// std::cout << "------------------------------------------------"  << endl;
		InputFile.close();
		
		OutputFile << "------------------------------------------------" << endl;
		// std::cout << "------------------------------------------------" << endl;
		OutputFile.close();
		
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
	}

	catch (TimeOutException& ) {
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
		if (extended_COV) {
			cov->add(cell->box);
		} else {
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
	// 	cout << " [total=" << cov->time() << "]";
	// cout << endl;
	// cout << " number of cells:\t\t" << nb_cells;
	// if (cov->nb_cells()!=nb_cells)
	// 	cout << " [total=" << cov->nb_cells() << "]";
	// cout << endl << endl;

	// if (statistics)
	// 	cout << "  ===== Statistics ====" << endl << endl << *statistics << endl;
	cout << nb_cells << " " << time << endl;
}



} // end namespace ibex